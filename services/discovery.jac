node Target {
    has name: str = "";
    has disease_context: str = "";
}

node Candidate {
    has smiles: str = "";
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
        candidate = Candidate(smiles=smiles, generation=generation, mol_weight=weight, logp=logp_value, h_donors=donor_count, h_acceptors=acceptor_count, lipinski_pass=passed, proxy_score=score, failure_reason=reason, mutation_applied=mutation);
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
        views.append(CandidateView(smiles=candidate.smiles, generation=candidate.generation, mol_weight=candidate.mol_weight, logp=candidate.logp, h_donors=candidate.h_donors, h_acceptors=candidate.h_acceptors, lipinski_pass=candidate.lipinski_pass, proxy_score=candidate.proxy_score, failure_reason=candidate.failure_reason, mutation_applied=candidate.mutation_applied));
    }
    return DiscoveryResult(target_name=target.name, disease_context=target.disease_context, candidates=views, docking_status="Not run", scope_note="This demo uses character-count descriptor estimates and scripted edit hypotheses. It does not resolve UniProt/PDB, validate SMILES with RDKit, use an LLM, or run AutoDock Vina.");
}
