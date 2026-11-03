def run(hops):
    return {
        "report": f"{OUT}/{{genome}}/{{assembly_type}}/eris/h{hops}_report.tsv",
        "cds": f"{OUT}/{{genome}}/{{assembly_type}}/eris/h{hops}_cds.fna",
    }


SCORED, REACH = run(config["hops"]), run(config["reach_hops"])


rule score:
    input:
        report=SCORED["report"],
        cds=SCORED["cds"],
        reach_report=REACH["report"],
        reach_cds=REACH["cds"],
        reference=lambda wc: sample(wc).reference,
        truth=lambda wc: sample(wc).truth,
        groups=rules.clone_eris.output.groups,
    output:
        multiext(OUT + "/{genome}/{assembly_type}/",
                 "metrics.tsv", "details.tsv", "relations.tsv", "side_confusion.tsv", "distance_confusion.tsv",
                 "effect_confusion.tsv"),
    params:
        outdir=OUT + "/{genome}/{assembly_type}",
        replicons=lambda wc: f"--replicons {sample(wc).replicons}" if sample(wc).replicons else "",
        thresholds=" ".join(f"--{k.replace('_', '-')} {config[k]}"
                            for k in ("min_overlap", "near", "min_identity", "min_cds_coverage")),
    log:
        OUT + "/{genome}/{assembly_type}/score.log",
    shell:
        "erisval score --report {input.report} --cds {input.cds}"
        " --reach-report {input.reach_report} --reach-cds {input.reach_cds}"
        " --reference {input.reference} --truth {input.truth} --groups {input.groups} {params.replicons}"
        " {params.thresholds} -o {params.outdir}"
        " --label genome={wildcards.genome} --label assembly_type={wildcards.assembly_type} 2> {log}"


rule combine:
    input:
        lambda wc: [f"{OUT}/{g}/{a}/{wc.table}.tsv" for g, a in SAMPLES.index],
    output:
        OUT + "/{table,metrics|details}.tsv",
    shell:
        "erisval combine -o {output} {input}"
