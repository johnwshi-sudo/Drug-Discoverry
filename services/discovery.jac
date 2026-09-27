node Target {
    has name: str = "";
    has disease_context: str = "";
}

node Candidate {
    has smiles: str = "";
    has molecular_formula: str = "";
    has generation: int = 0;
    has mol_weight: float = 0.0;
    has logp: float = 0.0;
    has h_donors: int = 0;
    has h_acceptors: int = 0;
    has lipinski_pass: bool = False;
    has proxy_score: float = 0.0;
    has failure_reason: str = "";
    has mutation_applied: str = "Seed molecule";
}

edge DerivedFrom: Candidate --> Candidate {
    has mutation_applied: str = "";
}

obj CandidateView {
    has smiles: str;
    has molecular_formula: str;
    has generation: int;
    has mol_weight: float;
    has logp: float;
    has h_donors: int;
    has h_acceptors: int;
    has lipinski_pass: bool;
    has proxy_score: float;
    has failure_reason: str;
    has mutation_applied: str;
}

obj DiscoveryResult {
    has target_name: str;
    has disease_context: str;
    has candidates: list[CandidateView];
    has docking_status: str;
    has scope_note: str;
}

def _estimate_weight(smiles: str) -> float {
    return float(len(smiles) * 12 + 96);
}

def _estimate_logp(smiles: str) -> float {
    return float(smiles.count("C") * 0.35 + smiles.count("c") * 0.18 - smiles.count("O") * 0.45 - smiles.count("N") * 0.25);
}

def _subscript_digits(value: int) -> str {
    result = "";
    for digit in str(value) {
        if digit == "0" { result += "₀"; }
        elif digit == "1" { result += "₁"; }
        elif digit == "2" { result += "₂"; }
        elif digit == "3" { result += "₃"; }
        elif digit == "4" { result += "₄"; }
        elif digit == "5" { result += "₅"; }
        elif digit == "6" { result += "₆"; }
        elif digit == "7" { result += "₇"; }
        elif digit == "8" { result += "₈"; }
        elif digit == "9" { result += "₉"; }
    }
    return result;
}

def _molecular_formula(smiles: str) -> str {
    # Parse the neutral organic subset needed by this demo. Unsupported or
    # malformed SMILES deliberately return an empty string instead of guessing.
    elements: list[str] = [];
    aromatic_atoms: list[bool] = [];
    bracketed_atoms: list[bool] = [];
    explicit_hydrogens: list[int] = [];
    bond_totals: list[float] = [];
    branches: list[int] = [];
    ring_atoms: dict[str, int] = {};
    ring_orders: dict[str, float] = {};
    current = -1;
    pending_order = 0.0;
    i = 0;
    element: str = "";

    while i < len(smiles) {
        char = smiles[i];
        if char == "[" {
            end = i + 1;
            while end < len(smiles) and smiles[end] != "]" { end += 1; }
            if end >= len(smiles) { return ""; }
            body = smiles[i + 1:end];
            if len(body) == 0 { return ""; }

            body_index = 0;
            if body[body_index].isdigit() { return ""; }
            atom_char = body[body_index];
            is_aromatic = atom_char.islower();
            element = atom_char.upper();
            body_index += 1;
            if not is_aromatic and body_index < len(body) and body[body_index].islower() {
                element += body[body_index];
                body_index += 1;
            }
            if element not in ["B", "C", "N", "O", "P", "S", "F", "Cl", "Br", "I"] {
                return "";
            }

            explicit_h = 0;
            while body_index < len(body) {
                part = body[body_index];
                if part == "@" {
                    while body_index < len(body) and body[body_index] == "@" { body_index += 1; }
                } elif part == "H" {
                    body_index += 1;
                    hydrogen_text = "";
                    while body_index < len(body) and body[body_index].isdigit() {
                        hydrogen_text += body[body_index];
                        body_index += 1;
                    }
                    explicit_h = 1 if hydrogen_text == "" else int(hydrogen_text);
                } else {
                    # Charge, isotope, atom-map, and other advanced bracket
                    # annotations require a chemistry toolkit for safe parsing.
                    return "";
                }
            }

            elements.append(element);
            aromatic_atoms.append(is_aromatic);
            bracketed_atoms.append(True);
            explicit_hydrogens.append(explicit_h);
            bond_totals.append(0.0);
            i = end + 1;
        } elif char == "(" {
            if current < 0 { return ""; }
            branches.append(current);
            i += 1;
            continue;
        } elif char == ")" {
            if len(branches) == 0 { return ""; }
            current = branches.pop();
            i += 1;
            continue;
        } elif char == "." {
            current = -1;
            i += 1;
            continue;
        } elif char in ["-", "=", "#", ":", "/", "\\"] {
            if current < 0 or pending_order != 0.0 { return ""; }
            if char == "=" { pending_order = 2.0; }
            elif char == "#" { pending_order = 3.0; }
            elif char == ":" { pending_order = 1.5; }
            else { pending_order = 1.0; }
            i += 1;
            continue;
        } elif char.isdigit() or char == "%" {
            if current < 0 { return ""; }
            label = char;
            if char == "%" {
                if i + 2 >= len(smiles) or not smiles[i + 1].isdigit() or not smiles[i + 2].isdigit() { return ""; }
                label = smiles[i + 1] + smiles[i + 2];
                i += 3;
            } else {
                i += 1;
            }
            prior_atom = ring_atoms.get(label, -1);
            if prior_atom < 0 {
                ring_atoms[label] = current;
                ring_orders[label] = pending_order;
                pending_order = 0.0;
                continue;
            }
            bond_order = pending_order;
            if bond_order == 0.0 { bond_order = ring_orders.get(label, 0.0); }
            if bond_order == 0.0 {
                bond_order = 1.5 if aromatic_atoms[current] and aromatic_atoms[prior_atom] else 1.0;
            }
            bond_totals[current] += bond_order;
            bond_totals[prior_atom] += bond_order;
            ring_atoms[label] = -1;
            ring_orders[label] = 0.0;
            pending_order = 0.0;
            continue;
        } else {
            element = "";
            is_aromatic = char.islower();
            if char in ["b", "c", "n", "o", "p", "s"] {
                element = char.upper();
            } elif char.isupper() {
                element = char;
                if i + 1 < len(smiles) and (char + smiles[i + 1]) in ["Cl", "Br"] {
                    element += smiles[i + 1];
                    i += 1;
                }
            } else {
                return "";
            }
            if element not in ["B", "C", "N", "O", "P", "S", "F", "Cl", "Br", "I"] { return ""; }
            explicit_h = 0;
            elements.append(element);
            aromatic_atoms.append(is_aromatic);
            bracketed_atoms.append(False);
            explicit_hydrogens.append(explicit_h);
            bond_totals.append(0.0);
            i += 1;
        }

        new_atom = len(elements) - 1;
        if current >= 0 {
            bond_order = pending_order;
            if bond_order == 0.0 {
                bond_order = 1.5 if aromatic_atoms[current] and aromatic_atoms[new_atom] else 1.0;
            }
            bond_totals[current] += bond_order;
            bond_totals[new_atom] += bond_order;
        }
        current = new_atom;
        pending_order = 0.0;
    }

    if len(elements) == 0 or len(branches) != 0 or pending_order != 0.0 { return ""; }
    for label in ring_atoms.keys() {
        if ring_atoms[label] >= 0 { return ""; }
    }

    counts: dict[str, int] = {};
    hydrogen_total = 0;
    for atom_index in range(len(elements)) {
        element = elements[atom_index];
        counts[element] = counts.get(element, 0) + 1;
        hydrogen_total += explicit_hydrogens[atom_index];
        if bracketed_atoms[atom_index] { continue; }

        valence = 0.0;
        if aromatic_atoms[atom_index] {
            valence = 4.0 if element == "C" else 3.0;
        } elif element == "C" or element == "Si" { valence = 4.0; }
        elif element == "B" { valence = 3.0; }
        elif element == "N" { valence = 3.0; }
        elif element == "O" or element == "Se" { valence = 2.0; }
        elif element == "P" { valence = 3.0 if bond_totals[atom_index] <= 3.0 else 5.0; }
        elif element == "S" { valence = 2.0 if bond_totals[atom_index] <= 2.0 else (4.0 if bond_totals[atom_index] <= 4.0 else 6.0); }
        else { valence = 1.0; }

        remaining = valence - bond_totals[atom_index];
        if remaining < 0.0 or float(int(remaining)) != remaining { return ""; }
        hydrogen_total += int(remaining);
    }
    if hydrogen_total > 0 { counts["H"] = hydrogen_total; }

    formula_order = ["B", "Br", "Cl", "F", "H", "I", "N", "O", "P", "S"];
    if counts.get("C", 0) > 0 {
        formula_order = ["C", "H", "B", "Br", "Cl", "F", "I", "N", "O", "P", "S"];
    }
    formula = "";
    for element in formula_order {
        count = counts.get(element, 0);
        if count > 0 {
            formula += element;
            if count > 1 { formula += _subscript_digits(count); }
        }
    }
    return formula;
}

def:pub run_discovery(target_name: str, disease_context: str, seed_smiles: str) -> DiscoveryResult {
    if not target_name.strip() or not seed_smiles.strip() {
        raise ValueError("Enter a target and a seed SMILES string.");
    }
    target = Target(name=target_name.strip(), disease_context=disease_context.strip());
    root ++> target;
    seed = seed_smiles.strip();
    generation_count = 4;
    candidates: list[Candidate] = [];
    previous: Candidate | None = None;
    for generation in range(generation_count) {
        smiles = seed;
        mutation = "Seed molecule";
        if generation == 1 {
            mutation = "Heuristic edit: hydroxyl substitution hypothesis";
            if "O" not in seed { smiles = seed + "O"; }
        } elif generation == 2 {
            mutation = "Heuristic edit: nitrogen substitution hypothesis";
            if "N" not in seed { smiles = seed + "N"; }
        } elif generation == 3 {
            mutation = "Heuristic edit: polarity adjustment hypothesis";
            smiles = seed + "O";
            if "N" not in smiles { smiles = smiles + "N"; }
        }
        weight = _estimate_weight(smiles);
        logp_value = _estimate_logp(smiles);
        donor_count = smiles.count("O") + smiles.count("N");
        acceptor_count = donor_count;
        passed = weight <= 500.0 and logp_value <= 5.0 and donor_count <= 5 and acceptor_count <= 10;
        score = 0.82;
        if weight > 500.0 { score = score - 0.18; }
        if logp_value > 5.0 { score = score - 0.16; }
        if donor_count > 5 { score = score - 0.12; }
        if acceptor_count > 10 { score = score - 0.10; }
        if score < 0.0 { score = 0.0; }
        reason = "Within heuristic Lipinski thresholds";
        if not passed { reason = "Heuristic descriptor estimate exceeds at least one Rule of Five threshold"; }
        candidate = Candidate(smiles=smiles, molecular_formula=_molecular_formula(smiles), generation=generation, mol_weight=weight, logp=logp_value, h_donors=donor_count, h_acceptors=acceptor_count, lipinski_pass=passed, proxy_score=score, failure_reason=reason, mutation_applied=mutation);
        root ++> candidate;
        if previous is not None {
            prior = previous as Candidate;
            prior +>:DerivedFrom(mutation_applied=mutation):+> candidate;
        }
        candidates.append(candidate);
        previous = candidate;
    }
    views: list[CandidateView] = [];
    for candidate in candidates {
        views.append(CandidateView(smiles=candidate.smiles, molecular_formula=candidate.molecular_formula, generation=candidate.generation, mol_weight=candidate.mol_weight, logp=candidate.logp, h_donors=candidate.h_donors, h_acceptors=candidate.h_acceptors, lipinski_pass=candidate.lipinski_pass, proxy_score=candidate.proxy_score, failure_reason=candidate.failure_reason, mutation_applied=candidate.mutation_applied));
    }
    return DiscoveryResult(target_name=target.name, disease_context=target.disease_context, candidates=views, docking_status="Not run", scope_note="This demo uses scripted edit hypotheses and heuristic descriptors. Its basic neutral-organic SMILES parser displays formulas when supported; it does not validate molecules with RDKit, use an LLM, or run docking.");
}
