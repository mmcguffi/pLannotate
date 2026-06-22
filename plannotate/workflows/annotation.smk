"""Fan out pLannotate database searches as independent local jobs."""

DATABASE_INDICES = range(int(config["database_count"]))


rule all:
    input:
        expand("results/{database_index}.pkl", database_index=DATABASE_INDICES)


rule search_database:
    output:
        "results/{database_index}.pkl"
    threads:
        lambda wildcards: config["thread_allocations"][int(wildcards.database_index)]
    params:
        python=config["python_executable"],
        query=config["query_file"],
        yaml=config["yaml_file"],
        linear="--linear" if config["linear"] else ""
    shell:
        """
        {params.python:q} -m plannotate._annotation_worker \
            --query {params.query:q} \
            --yaml {params.yaml:q} \
            --database-index {wildcards.database_index} \
            --output {output:q} \
            --threads {threads} \
            {params.linear}
        """
