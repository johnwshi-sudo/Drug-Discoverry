sv import from ..services.discovery { CandidateView, run_discovery }
import from .ui.button { Button }
import from .ui.card { Card, CardContent, CardDescription, CardHeader, CardTitle }
import from .ui.input { Input }
import from .ui.field { Field, FieldLabel }
import from .ui.badge { Badge }
cl import from "@hugeicons/react" { HugeiconsIcon }
cl import from "@hugeicons/core-free-icons" { SparklesIcon, Atom01Icon, Alert01Icon, CheckmarkCircle01Icon, ArrowDown01Icon }

def:pub CandidateChain(candidates: list[CandidateView]) -> JsxElement {
    return <Card className="shadow-sm"><CardHeader><CardTitle>Refinement chain</CardTitle><CardDescription>Graph-linked candidate hypotheses · heuristic descriptor estimates</CardDescription></CardHeader><CardContent className="flex flex-col gap-4">
        {for (index, candidate) in enumerate(candidates) {
            <article key={str(index)} className="relative rounded-xl border bg-background p-4 shadow-sm">
                {if index < len(candidates) - 1 {<div className="absolute -bottom-5 left-7 z-10 flex size-6 items-center justify-center rounded-full border bg-background text-muted-foreground"><HugeiconsIcon icon={ArrowDown01Icon} strokeWidth={2} className="size-3" /></div>}}
                <div className="flex flex-wrap items-center justify-between gap-2"><div className="flex items-center gap-2"><span className="flex size-7 items-center justify-center rounded-full bg-primary/10 text-xs font-semibold text-primary">{str(candidate.generation)}</span><h3 className="text-sm font-semibold">Generation {str(candidate.generation)}</h3></div><Badge variant={"secondary" if candidate.lipinski_pass else "outline"}>{"Proxy score " + str(candidate.proxy_score)}</Badge></div>
                <code className="mt-3 block break-all rounded-md bg-muted/70 p-3 font-mono text-xs leading-relaxed">{candidate.smiles}</code><p className="mt-3 text-xs text-muted-foreground">{candidate.mutation_applied}</p>
                <div className="mt-3 grid grid-cols-2 gap-2 sm:grid-cols-4"><div className="rounded-lg bg-muted/50 p-2"><p className="text-[11px] text-muted-foreground">MW estimate</p><p className="mt-1 text-sm font-medium tabular-nums">{str(candidate.mol_weight)}</p></div><div className="rounded-lg bg-muted/50 p-2"><p className="text-[11px] text-muted-foreground">LogP estimate</p><p className="mt-1 text-sm font-medium tabular-nums">{str(candidate.logp)}</p></div><div className="rounded-lg bg-muted/50 p-2"><p className="text-[11px] text-muted-foreground">H-bond donors</p><p className="mt-1 text-sm font-medium tabular-nums">{str(candidate.h_donors)}</p></div><div className="rounded-lg bg-muted/50 p-2"><p className="text-[11px] text-muted-foreground">H-bond acceptors</p><p className="mt-1 text-sm font-medium tabular-nums">{str(candidate.h_acceptors)}</p></div></div>
                <p className="mt-3 flex items-center gap-1.5 text-xs text-muted-foreground"><HugeiconsIcon icon={CheckmarkCircle01Icon if candidate.lipinski_pass else Alert01Icon} strokeWidth={2} className="size-3.5" />{candidate.failure_reason}</p>
            </article>
        }}
    </CardContent></Card>;
}

def:pub DiscoveryWorkspace() -> JsxElement {
    has target_name: str = "EGFR kinase";
    has disease_context: str = "Non-small cell lung cancer";
    has seed_smiles: str = "CCOc1ccc(NC(=O)C)cc1";
    has candidates: list[CandidateView] = [];
    has run_target: str = "";
    has docking_status: str = "";
    has scope_note: str = "";
    has loading: bool = False;
    has error_message: str = "";

    def change_target(e: ChangeEvent) { target_name = e.target.value; }
    def change_context(e: ChangeEvent) { disease_context = e.target.value; }
    def change_smiles(e: ChangeEvent) { seed_smiles = e.target.value; }
    async def handle_run -> None {
        loading = True;
        error_message = "";
        try {
            outcome = await run_discovery(target_name, disease_context, seed_smiles);
            candidates = outcome.candidates;
            run_target = outcome.target_name + " · " + outcome.disease_context;
            docking_status = outcome.docking_status;
            scope_note = outcome.scope_note;
        } except Exception {
            error_message = "The run could not be completed. Check the target and SMILES input, then try again.";
        }
        loading = False;
    }

    return <main className="min-h-screen bg-muted/40"><header className="border-b bg-background"><div className="mx-auto flex max-w-7xl items-center justify-between gap-4 px-5 py-4 sm:px-8"><div className="flex items-center gap-3"><div className="flex size-10 items-center justify-center rounded-xl bg-primary text-primary-foreground"><HugeiconsIcon icon={Atom01Icon} strokeWidth={2} className="size-5" /></div><div><p className="font-semibold tracking-tight">Asterion Discovery</p><p className="text-xs text-muted-foreground">Molecular research workspace</p></div></div><Badge variant="outline">Proof of concept</Badge></div></header>
        <div className="mx-auto flex w-full max-w-7xl flex-col gap-6 px-5 py-8 sm:px-8"><div className="flex flex-wrap items-end justify-between gap-3"><div className="flex flex-col gap-1"><h1 className="text-2xl font-semibold tracking-tight">Candidate refinement</h1><p className="text-sm text-muted-foreground">Explore a graph-backed, heuristic refinement demo against a molecular target.</p></div><span className="text-xs text-muted-foreground">Iterative hypothesis chain · up to 6 generations</span></div>
            <div className="grid gap-6 lg:grid-cols-[minmax(300px,0.78fr)_minmax(0,1.55fr)]"><Card className="h-fit shadow-sm"><CardHeader><CardTitle>Research setup</CardTitle><CardDescription>Choose a target context and provide a seed molecule.</CardDescription></CardHeader><CardContent className="flex flex-col gap-5"><Field><FieldLabel htmlFor="target-name">Target name</FieldLabel><Input id="target-name" value={target_name} onChange={change_target} placeholder="EGFR kinase" /></Field><Field><FieldLabel htmlFor="disease-context">Disease context</FieldLabel><Input id="disease-context" value={disease_context} onChange={change_context} placeholder="Disease or biological context" /></Field><Field><FieldLabel htmlFor="seed-smiles">Seed SMILES</FieldLabel><textarea id="seed-smiles" value={seed_smiles} onChange={change_smiles} placeholder="Enter a molecular SMILES string" className="min-h-24 w-full rounded-md border border-input bg-background px-3 py-2 font-mono text-sm shadow-xs outline-none focus-visible:ring-2 focus-visible:ring-ring" /></Field><p className="text-xs leading-relaxed text-muted-foreground">Demo input only. SMILES validity is not checked by a chemistry toolkit.</p><Button className="w-full" disabled={loading} onClick={handle_run}><HugeiconsIcon icon={SparklesIcon} strokeWidth={2} className="size-4" />{("Running refinement…" if loading else "Run refinement demo")}</Button>{error_message and <p role="alert" className="text-sm text-destructive">{error_message}</p>}</CardContent></Card>
            <section className="flex min-w-0 flex-col gap-5">{if len(candidates) == 0 {<Card className="shadow-sm"><CardContent className="flex min-h-72 flex-col items-center justify-center gap-3 text-center"><div className="flex size-12 items-center justify-center rounded-xl bg-primary/10 text-primary"><HugeiconsIcon icon={Atom01Icon} strokeWidth={2} className="size-6" /></div><h2 className="text-lg font-semibold">Your candidate chain starts here</h2><p className="max-w-md text-sm text-muted-foreground">Set a target and seed molecule, then run the demo to inspect each graph-linked hypothesis and its estimated descriptors.</p></CardContent></Card>}}
                {if len(candidates) > 0 {<><div className="flex flex-wrap items-center justify-between gap-2"><h2 className="text-sm font-medium">{run_target}</h2><Badge variant="outline">Heuristic proxy</Badge></div><CandidateChain candidates={candidates} /><Card className="border-dashed shadow-sm"><CardHeader><CardTitle className="text-base">Docking status: {docking_status}</CardTitle><CardDescription>{scope_note}</CardDescription></CardHeader></Card></>}}</section></div>
            <p className="text-center text-xs leading-relaxed text-muted-foreground">Research prototype only. Heuristic descriptor estimates are not RDKit results, clinical evidence, or validated drug-discovery recommendations.</p></div></main>;
}
