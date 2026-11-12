const fs = require("fs");
const d = require("docx");
const { Document, Packer, Paragraph, TextRun, HeadingLevel, AlignmentType,
        Table, TableRow, TableCell, WidthType, ShadingType, BorderStyle,
        ImageRun, PageBreak } = d;

const F = "Times New Roman";
const INK = "000000", MUTE = "444444", ACC = "1F3864";

const P = (t, o = {}) => new Paragraph({
  spacing: { after: o.after === undefined ? 140 : o.after, line: 320 },
  alignment: o.align || AlignmentType.JUSTIFIED, indent: o.indent,
  children: [new TextRun({ text: t, size: o.size || 24, font: F,
    color: o.color || INK, bold: o.bold, italics: o.italics })]
});
const PR = (runs, o = {}) => new Paragraph({
  spacing: { after: o.after === undefined ? 140 : o.after, line: 320 },
  alignment: o.align || AlignmentType.JUSTIFIED,
  children: runs.map(r => new TextRun({ text: r.t, bold: r.b, italics: r.i,
    color: r.c || INK, size: o.size || 24, font: F }))
});
const H1 = t => new Paragraph({ text: t, heading: HeadingLevel.HEADING_1 });
const H2 = t => new Paragraph({ text: t, heading: HeadingLevel.HEADING_2 });
const CODE = lines => new Paragraph({
  spacing: { before: 80, after: 140, line: 240 },
  shading: { type: ShadingType.CLEAR, fill: "F4F4F4", color: "auto" },
  border: { left: { style: BorderStyle.SINGLE, size: 12, color: "B8C6D4", space: 6 } },
  indent: { left: 240 },
  children: lines.flatMap((l, i) => [
    ...(i ? [new TextRun({ break: 1 })] : []),
    new TextRun({ text: l, font: "Consolas", size: 19, color: "1a1a1a" })])
});
const cell = (t, o = {}) => new TableCell({
  width: { size: o.w, type: WidthType.DXA },
  shading: o.shade ? { type: ShadingType.CLEAR, fill: o.shade, color: "auto" } : undefined,
  margins: { top: 70, bottom: 70, left: 110, right: 110 },
  verticalAlign: d.VerticalAlign.CENTER,
  children: (Array.isArray(t) ? t : [t]).map(x => new Paragraph({
    spacing: { after: 0, line: 260 }, alignment: o.align || AlignmentType.LEFT,
    children: [new TextRun({ text: x, bold: o.b, italics: o.i, size: o.size || 20,
      font: F, color: o.color || INK })] }))
});
function table(widths, header, rows) {
  return new Table({
    columnWidths: widths,
    width: { size: widths.reduce((a, b) => a + b, 0), type: WidthType.DXA },
    borders: {
      top: { style: BorderStyle.SINGLE, size: 4, color: "999999" },
      bottom: { style: BorderStyle.SINGLE, size: 4, color: "999999" },
      left: { style: BorderStyle.SINGLE, size: 4, color: "999999" },
      right: { style: BorderStyle.SINGLE, size: 4, color: "999999" },
      insideHorizontal: { style: BorderStyle.SINGLE, size: 2, color: "CCCCCC" },
      insideVertical: { style: BorderStyle.SINGLE, size: 2, color: "CCCCCC" } },
    rows: [
      new TableRow({ tableHeader: true, children: header.map((h, i) =>
        cell(h, { w: widths[i], b: true, shade: "1F3864", color: "FFFFFF", size: 20 })) }),
      ...rows.map((r, ri) => new TableRow({ cantSplit: true,
        children: r.map((c, i) => cell(c, { w: widths[i],
          shade: ri % 2 ? "F2F4F8" : undefined, size: 20 })) }))]
  });
}
const CAP = t => new Paragraph({
  spacing: { before: 110, after: 240 }, alignment: AlignmentType.CENTER,
  children: [new TextRun({ text: t, size: 20, italics: true, font: F, color: MUTE })]
});

const img = fs.readFileSync("workflow.png");
const B = [];

B.push(new Paragraph({ spacing: { before: 420, after: 80 }, alignment: AlignmentType.CENTER,
  children: [new TextRun({ text: "Construction of the Anticancer Natural Product Screening Library",
    bold: true, size: 34, font: F, color: ACC })] }));
B.push(new Paragraph({ spacing: { after: 320 }, alignment: AlignmentType.CENTER,
  children: [new TextRun({ text: "Source database, procedure and rationale for virtual screening against PARP7",
    size: 25, font: F, color: MUTE })] }));

B.push(H1("1.  Scope"));
B.push(P("This document sets out how the library of anticancer natural products used in this work was assembled for structure based virtual screening against the PARP7 catalytic domain model. It describes the source database, gives the procedure as a sequence of executable steps, states the reasoning behind every narrowing decision, and explains why a library of this kind is well suited to this particular target. The accompanying code repository reproduces the procedure in full."));

B.push(H1("2.  The source database"));
B.push(P("The library was assembled from NPASS, the Natural Product Activity and Species Source database (Zhao et al., 2023). The release used holds 204,023 natural products, 1,048,756 quantitative activity records and 8,764 biological targets, of which 1,658 are human cancer cell lines. Structures, activity records and target descriptions are distributed as separate tables, so the anticancer subset was assembled by joining them."));
B.push(P("What distinguishes NPASS from a general natural product collection is that the anticancer evidence it carries is quantitative and referenced rather than categorical. Every retained compound has a measured potency value, the cell line against which that value was obtained, and a literature citation. The consequence for this work is that the operational definition of anticancer activity is set by the investigator and can be stated in a single sentence, rather than being inherited from the inclusion criteria of a curator whose reasoning is not recorded."));
B.push(P("Of the 1,048,756 activity records in the release, 279,591 were measured against cancer cell lines. The potency threshold and the curation described below reduce these to the final library."));

B.push(H1("3.  Why this library suits PARP7"));
B.push(H2("3.1  The explored chemical space is unusually narrow"));
B.push(P("Almost every reported PARP7 inhibitor descends from a single clinical compound, RBN-2397, and its close analogues, which are nicotinamide mimetic heterocycles (Gozgit et al., 2021). For a target whose entire body of chemical matter comes from one synthetic lineage, the most useful contribution a virtual screen can make is scaffold novelty. Natural products supply ring systems, substitution patterns and stereochemical arrangements that this medicinal chemistry programme has never visited."));
B.push(H2("3.2  The binding site recognises a natural metabolite"));
B.push(P("The PARP7 catalytic site binds NAD+, and its nicotinamide subpocket recognises a small endogenous metabolite. Natural products are biosynthesised by enzymes that handle these same cofactors, and heteroaromatic amide motifs resembling nicotinamide are common among them. A population of molecules shaped by selection to interact with NAD+ handling machinery is a more plausible starting point for this pocket than a synthetic screening collection assembled without reference to it."));
B.push(H2("3.3  Two independent lines of evidence converge on each hit"));
B.push(P("A compound that both scores well against the receptor model and already carries a measured cancer cell line potency is a stronger candidate for follow up than one supported by a docking score alone. The compound is obtainable, and the recorded potency gives a concentration range for assay design, which shortens the path from a computational hit to an experiment considerably."));
B.push(H2("3.4  The screen can raise a testable mechanistic hypothesis"));
B.push(P("PARP7 inhibition acts by restoring type I interferon signalling rather than by direct cytotoxicity (Gozgit et al., 2021). Many compounds in this library inhibit the growth of cancer cells by mechanisms that were never determined. Where such a compound docks well at the PARP7 NAD+ site, the result is a specific and testable proposal about a mechanism that no one has assigned. This is a question a phenotype annotated library can raise and a mechanism agnostic one cannot."));

B.push(new Paragraph({ children: [new PageBreak()] }));
B.push(H1("4.  Workflow"));
B.push(P("Figure 1 shows the procedure from the source database to the prepared ligands. Each step is described in Section 5 and the curation applied at step 4 is set out in Section 6."));
B.push(new Paragraph({ spacing: { before: 140, after: 0 }, alignment: AlignmentType.CENTER,
  children: [new ImageRun({ data: img, type: "png", transformation: { width: 520, height: 608 } })] }));
B.push(CAP("Figure 1.  Construction of the anticancer natural product screening library. The benchmark actives used to validate the docking protocol are drawn separately and are not part of this library, for the reason given in Section 7."));

B.push(H1("5.  Procedure"));

B.push(H2("Step 0.  Environment check"));
B.push(P("The environment was verified before any data was processed. The check confirms every required package and then runs the whole pipeline on built in molecules, asserting that a malformed structure is rejected, that a sugar only molecule is removed, that rutin is deglycosylated and merges onto quercetin, that quercetin is retained but flagged, that a lipophilic triterpene survives the drug likeness filter, that a counterion is stripped, and that the prepared ligands carry genuine three dimensional coordinates and a torsion tree. An environment problem is therefore found in seconds rather than part way through a long run."));
B.push(CODE(["conda env create -f environment.yml && conda activate deck",
             "python scripts/00_preflight.py"]));

B.push(H2("Step 1.  Acquisition"));
B.push(P("The NPASS structure, activity, target and general information tables were downloaded, and the release recorded. The download is verified rather than assumed, since a site that answers a missing file with an error page still returns a success code and the resulting file is HTML rather than data."));
B.push(CODE(["bash scripts/00_fetch_sources.sh"]));

B.push(H2("Step 2.  Selection of cancer cell-line actives"));
B.push(P("The activity table was joined to the target table and the records whose target is a human cancer cell line were retained. Activity values were normalised to micromolar before the threshold was applied, so that records stored in nanomolar or molar are compared on the same scale, and where a compound carries several records the most potent was kept. Records reported in mass units such as micrograms per millilitre cannot be converted without a molecular weight and were retained with a flag rather than discarded."));
B.push(CODE(["python scripts/11_npass_join.py \\",
             "    --activities data/raw/npass_raw/npass_activities.txt \\",
             "    --structures data/raw/npass_raw/npass_structure.txt \\",
             "    --targets    data/raw/npass_raw/npass_target.txt \\",
             "    --max-um 10 \\",
             "    --out data/raw/npass_anticancer.csv"]));
B.push(PR([{ t: "The threshold of 10 micromolar is the operational definition of anticancer activity used in this work. ", b: true },
  { t: "It is a conventional cutoff for a cell line screen and is recorded here so that the definition is explicit rather than implied. The activity types accepted were IC50, GI50, EC50, CC50, LC50 and ED50." }]));

B.push(H2("Step 3.  Assembly and merging"));
B.push(P("The selected compounds were assembled into a single file. Where further anticancer compound lists are available they are read at this point and merged, with the column holding the structure detected automatically. Compounds are deduplicated across sources on full InChIKey and each records the sources in which it occurred, so that a compound seen in more than one database can be distinguished from one seen once."));
B.push(CODE(["python scripts/10_build_from_sources.py --auto data/raw \\",
             "    --out data/anticancer_raw.sdf --report logs/source_overlap.tsv"]));

B.push(H2("Step 4.  Chemical curation"));
B.push(P("The curation set out in Section 6 was applied. Every compound removed was written to a rejection log together with the reason."));
B.push(CODE(["python scripts/01_filter.py --in data/anticancer_raw.sdf \\",
             "    --out data/anticancer_final.sdf \\",
             "    --props logs/props.tsv --rejects logs/rejects.tsv"]));

B.push(H2("Step 5.  Three dimensional preparation"));
B.push(P("Three dimensional coordinates were generated once, at this stage and nowhere earlier, since nothing in the curation requires them. Hydrogens were added at pH 7.4, coordinates embedded with ETKDGv3, geometry optimised with MMFF94, and each ligand written to PDBQT with an explicit torsion tree by Meeko."));
B.push(P("Protonation was applied here rather than during standardisation because the relevant protonation state is a property of the assay conditions and not of the compound. Meeko was used in preference to a generic format converter because it writes the torsion tree that AutoDock reads to determine which bonds rotate. A file converted without a correct torsion tree will dock, but with the wrong degrees of freedom, and the failure is silent."));
B.push(CODE(["python scripts/04_prepare_pdbqt.py --in data/anticancer_final.sdf \\",
             "    --out data/pdbqt --jobs 0 --failures logs/pdbqt_failures.tsv"]));

B.push(H2("Step 6.  Export and provenance"));
B.push(P("The assembled and curated compound sets were exported to SMILES for deposition, together with a table of the per compound annotations, the software versions, the checksum of the library and the attrition by reason."));
B.push(CODE(["bash run_all.sh      # runs steps 1 to 6 with a check after each"]));

B.push(new Paragraph({ children: [new PageBreak()] }));
B.push(H1("6.  Chemical curation"));
B.push(P("Each step is given with the criterion applied and the reason it is justified. Every compound removed at every step was written to a rejection log together with the reason, so that the attrition can be audited and a discarded subset recovered if a later question requires it."));

B.push(table([1750, 2500, 4776],
  ["Step", "Criterion", "Reason"],
  [
    ["Sanitisation", "Valence and aromaticity perceived successfully", "Every property and conformer computed later depends on a correct valence model. A record that cannot be sanitised has no meaningful molecular weight and no meaningful docking score."],
    ["Standardisation", "Largest fragment retained, counterions removed, charge neutralised", "A counterion is a property of the isolated sample rather than of the molecule that binds the target. Neutralisation also gives a single canonical form so that duplicate removal compares molecules rather than salt forms."],
    ["Element whitelist", "H, B, C, N, O, F, Si, P, S, Cl, Br and I only", "Empirical scoring functions carry no parameters for transition metal coordination. A score assigned to an organometallic compound is the output of a formula applied outside the domain in which it was fitted."],
    ["Deglycosylation", "Terminal circular sugars cleaved, sugar only molecules removed", "Glycosidic chains inflate the rotatable bond count and present hydroxyls that scoring functions reward spuriously. The aglycone is frequently the circulating species, since dietary glycosides are hydrolysed at the intestinal brush border, and the PARP7 pocket is too small for a bulky glycoside."],
    ["Molecular weight", "150 to 600 Da, applied after deglycosylation", "The upper bound is set by the target rather than by convention, the reference inhibitor occupying the site at 523 Da. Applying the window after deglycosylation allows a large saponin to be reduced into range rather than discarded."],
    ["Drug likeness", "At most one Lipinski violation", "One violation is tolerated specifically to retain lipophilic terpenoids. A stricter filter would preferentially delete the scaffold classes that characterise plant secondary metabolism."],
    ["PAINS", "Flagged and retained, not removed", "PAINS matches a large proportion of natural products, quercetin among them. Removing matches outright deletes much of the chemistry the screen exists to interrogate, and deletes it invisibly."],
    ["Duplicate removal", "Unique full InChIKey", "Performed after standardisation and deglycosylation, because both create duplicates that do not exist in the source. The full key is used so that stereoisomers remain distinct."]
  ]));
B.push(CAP("Table 1.  Curation steps, the criterion applied at each and the reason for it."));

B.push(H1("7.  What this library is not"));
B.push(PR([{ t: "The library is not a set of benchmark actives and was not used as one. ", b: true },
  { t: "A compound carrying a reported cancer cell line IC50 is known to inhibit the growth of cancer cells by an unspecified mechanism. The overwhelming majority of such compounds do not act through PARP7, since cytotoxicity in a cell line is a phenotypic observation compatible with hundreds of molecular mechanisms. Employing these compounds as the positive class in a retrospective enrichment calculation against the PARP7 model would produce a receiver operating characteristic curve, and an area beneath it, that measures nothing about the ability of the model to recognise PARP7 binders." }]));
B.push(P("Actives for retrospective validation were drawn instead from compounds with measured activity against PARP7 itself, obtained from ChEMBL and BindingDB for UniProt Q7Z3E1 and from the published inhibitor series, and were screened against property matched decoys. That benchmark was constructed separately."));
B.push(P("Four further limitations should be stated. The receptor is a comparative model built upon a mutated PARP12 template, since no experimental structure of PARP7 exists, so docking scores are conditional upon that model and absolute binding energies should not be reported as predictions of affinity. Deglycosylation restricts the screen to aglycones, so glycosides whose sugar is essential to activity lie outside its scope by construction. The PAINS flag identifies substructures associated with assay interference in synthetic libraries and is a prompt for scrutiny rather than a verdict on any individual natural product. Finally, where a diversity selection has been applied, chemical space has been sampled rather than enumerated, so a compound absent from the library has not been tested and no conclusion regarding it is warranted."));
B.push(P("One consideration specific to this target should also be recorded. The nicotinamide subpocket is strongly conserved across the PARP family, so hits arising from this screen were counter docked against PARP1, PARP2 and PARP12 before any claim of selectivity was made."));

B.push(H1("8.  Reproducing this work"));
B.push(P("The code repository accompanying this document reproduces the procedure in full. The environment is checked before any data is processed, and the check runs the whole pipeline on built in molecules and asserts that the chemistry is correct, so that an environment problem is found in seconds rather than part way through a long run."));
B.push(CODE(["conda env create -f environment.yml && conda activate deck",
             "python scripts/00_preflight.py        # must pass before proceeding",
             "bash scripts/00_fetch_sources.sh",
             "bash run_all.sh"]));
B.push(P("Two files are deposited as the compound set itself. The first holds every compound assembled before curation and the second holds the curated library, both in SMILES format with one compound per line and its identifier alongside. A third file carries the per compound annotations, so that the chemistry and its provenance can both be deposited without either file becoming unreadable."));

B.push(table([2600, 2400, 4026],
  ["File", "Content", "Purpose"],
  [
    ["total_molecules.smi", "All compounds assembled from the sources", "The starting set, before any curation"],
    ["final_molecules.smi", "The curated library", "The compound set actually screened"],
    ["final_molecules_annotations.tsv", "InChIKey, sources, activity, PAINS flag, sugars removed", "Provenance and annotation for every compound"],
    ["logs/provenance.txt", "Counts, software versions, checksum", "Identifies the exact build"],
    ["logs/attrition_summary.txt", "Compounds removed by each criterion", "The audit of the curation"],
    ["logs/rejects.tsv", "Every removed compound with its reason", "Allows any discarded subset to be recovered"]
  ]));
B.push(CAP("Table 2.  Files deposited in the repository."));

B.push(H1("9.  Software"));
B.push(table([2600, 1600, 4826],
  ["Software", "Version", "Use"],
  [
    ["RDKit", "2026.09.1", "Parsing, standardisation, descriptors, PAINS catalogue, fingerprints, diversity selection"],
    ["Meeko", "0.8.0", "PDBQT preparation with explicit torsion trees"],
    ["AutoDock Vina", "1.2 series", "Docking engine for which the library was prepared"],
    ["Python", "3.11", "Scripting environment"],
    ["NumPy, pandas, SciPy", "current", "Numerical analysis and table handling"]
  ]));
B.push(CAP("Table 3.  Software used."));

B.push(H1("10.  References"));
[
 "Baell JB, Holloway GA. New substructure filters for removal of pan assay interference compounds (PAINS) from screening libraries and for their exclusion in bioassays. Journal of Medicinal Chemistry 2010;53:2719-2740.",
 "Butina D. Unsupervised data base clustering based on Daylight's fingerprint and Tanimoto similarity. Journal of Chemical Information and Computer Sciences 1999;39:747-750.",
 "Choi H, Cho SY, Pak HJ, Kim Y, Choi J, Lee YJ, Gong BH, Kang YS, Han T, Choi G, Cho Y, Lee S, Ryoo D, Park H. NPCARE: database of natural products and fractional extracts for cancer regulation. Journal of Cheminformatics 2017;9:2.",
 "Day AJ, Canada FJ, Diaz JC, Kroon PA, Mclauchlan R, Faulds CB, Plumb GW, Morgan MR, Williamson G. Dietary flavonoid and isoflavone glycosides are hydrolysed by the lactase site of lactase phlorizin hydrolase. FEBS Letters 2000;468:166-170.",
 "Eberhardt J, Santos-Martins D, Tillack AF, Forli S. AutoDock Vina 1.2.0: new docking methods, expanded force field, and Python bindings. Journal of Chemical Information and Modeling 2021;61:3891-3898.",
 "Gozgit JM, Vasbinder MM, Abo RP, et al. PARP7 negatively regulates the type I interferon response in cancer cells and its inhibition triggers antitumor immunity. Cancer Cell 2021;39:1214-1226.e10.",
 "Heller SR, McNaught A, Pletnev I, Stein S, Tchekhovskoi D. InChI, the IUPAC International Chemical Identifier. Journal of Cheminformatics 2015;7:23.",
 "Lipinski CA, Lombardo F, Dominy BW, Feeney PJ. Experimental and computational approaches to estimate solubility and permeability in drug discovery and development settings. Advanced Drug Delivery Reviews 1997;23:3-25.",
 "Mangal M, Sagar P, Singh H, Raghava GPS, Agarwal SM. NPACT: Naturally Occurring Plant-based Anti-cancer Compound-Activity-Target database. Nucleic Acids Research 2013;41(D1):D1124-D1129.",
 "Nemeth K, Plumb GW, Berrin JG, Juge N, Jacob R, Naim HY, Williamson G, Swallow DM, Kroon PA. Deglycosylation by small intestinal epithelial cell beta-glucosidases is a critical step in the absorption and metabolism of dietary flavonoid glycosides in humans. European Journal of Nutrition 2003;42:29-42.",
 "Rogers D, Hahn M. Extended-connectivity fingerprints. Journal of Chemical Information and Modeling 2010;50:742-754.",
 "Schaub J, Zielesny A, Steinbeck C, Sorokina M. Too sweet: cheminformatics for deglycosylation in natural products. Journal of Cheminformatics 2020;12:67.",
 "Veber DF, Johnson SR, Cheng HY, Smith BR, Ward KW, Kopple KD. Molecular properties that influence the oral bioavailability of drug candidates. Journal of Medicinal Chemistry 2002;45:2615-2623.",
 "Zhao H, Yang Y, Wang S, et al. NPASS database update 2023: quantitative natural product activity and species source database for biomedical research. Nucleic Acids Research 2023;51(D1):D621-D628."
].forEach(r => B.push(new Paragraph({
  spacing: { after: 120, line: 300 }, indent: { left: 420, hanging: 420 },
  alignment: AlignmentType.JUSTIFIED,
  children: [new TextRun({ text: r, size: 21, font: F, color: INK })] })));

const doc = new Document({
  creator: "PARP7 virtual screening project",
  title: "Construction of the Anticancer Natural Product Screening Library",
  styles: { default: {
    document: { run: { font: F, size: 24, color: INK } },
    heading1: { run: { font: F, size: 30, bold: true, color: ACC },
                paragraph: { spacing: { before: 360, after: 160 }, keepNext: true } },
    heading2: { run: { font: F, size: 25, bold: true, color: "2E5F8A" },
                paragraph: { spacing: { before: 260, after: 110 }, keepNext: true } },
    heading3: { run: { font: F, size: 24, bold: true, color: "2E5F8A" },
                paragraph: { spacing: { before: 200, after: 90 }, keepNext: true } } } },
  sections: [{ properties: { page: { margin: { top: 1440, right: 1440, bottom: 1440, left: 1440 } } },
               children: B }]
});
Packer.toBuffer(doc).then(b => {
  fs.writeFileSync("Anticancer_screening_library_procedure.docx", b);
  console.log("wrote", b.length, "bytes");
});
