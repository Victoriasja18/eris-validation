rule clone_eris:
    output:
        commit=f"{ERIS_DIR}/COMMIT",
        groups=f"{ERIS_DIR}/src/tests/isoform_groups.tsv",
    params:
        repo=config["eris"]["repo"],
        src=f"{ERIS_DIR}/src",
    log:
        f"{ERIS_DIR}/clone.log",
    shell:
        "(rm -rf {params.src} && git clone --quiet {params.repo} {params.src}"
        " && git -C {params.src} checkout --quiet {COMMIT}"
        " && git -C {params.src} rev-parse HEAD > {output.commit}) > {log} 2>&1"


rule install_eris:
    input:
        rules.clone_eris.output.commit,
    output:
        f"{ERIS_DIR}/venv/bin/eris",
    params:
        src=f"{ERIS_DIR}/src",
        venv=f"{ERIS_DIR}/venv",
    log:
        f"{ERIS_DIR}/install.log",
    shell:
        "(python -m venv {params.venv} && {params.venv}/bin/pip install --quiet {params.src}) > {log} 2>&1"


rule check_is_db:
    input:
        config["is_db"],
    output:
        "resources/is_db.ok",
    shell:
        "erisval check-is-db --db {input} > {output}"


rule run_eris:
    input:
        eris=rules.install_eris.output,
        is_db_ok=rules.check_is_db.output,
        assembly=lambda wc: sample(wc).assembly,
    output:
        report="results/{genome}/{assembly_type}/eris/h{hops}_report.tsv",
        cds="results/{genome}/{assembly_type}/eris/h{hops}_cds.fna",
    params:
        prefix="results/{genome}/{assembly_type}/eris/h{hops}",
        is_db=config["is_db"],
        mode=config["mode"],
        min_coverage=config["min_coverage"],
    threads: config["threads"]
    log:
        "results/{genome}/{assembly_type}/eris/h{hops}.log",
    shell:
        "{input.eris} -i {input.assembly} -d {params.is_db} -o {params.prefix} --mode {params.mode}"
        " --hops {wildcards.hops} -c {params.min_coverage} -t {threads} > {log} 2>&1"
