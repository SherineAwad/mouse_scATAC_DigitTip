configfile: "config.yaml"

SAMPLES = config["samples"]

COMPARISONS = {
    comparison["name"]: comparison
    for comparison in config["comparisons"]
}

ANNOTATE = config["annotate"]
MOTIF = config["motif"]


rule all:
    input:
        "analysed.h5ad",

        # Differential peak results
        expand(
            "{comparison}_significant_results.csv",
            comparison=COMPARISONS.keys()
        ),

        expand(
            "{comparison}_all_results.csv",
            comparison=COMPARISONS.keys()
        ),

        # Annotated peaks
        expand(
            "{comparison}_annotated.csv",
            comparison=ANNOTATE
        ),

        # Motif enrichment
        expand(
            "motifs_{comparison}_gain_motif_enrichment.csv",
            comparison=MOTIF
        ),

        expand(
            "motifs_{comparison}_loss_motif_enrichment.csv",
            comparison=MOTIF
        )


rule preprocess:
    input:
        config["samples_file"]

    output:
        expand(
            "{sample}.h5ad",
            sample=SAMPLES
        )

    shell:
        """
        python src/preprocess.py \
            --samples {input}
        """


rule filter:
    input:
        "{sample}.h5ad"

    output:
        "{sample}_filtered.h5ad"

    shell:
        """
        python src/filter.py \
            --input {input} \
            --output {output} \
            --prefix {wildcards.sample}
        """


rule merge:
    input:
        expand(
            "{sample}_filtered.h5ad",
            sample=SAMPLES
        )

    output:
        "merged.h5ad"

    shell:
        """
        python src/merge.py \
            --samples {config[samples_file]} \
            --output {output}
        """


rule removeDoublets:
    input:
        "merged.h5ad"

    output:
        "noDoublets.h5ad"

    shell:
        """
        python src/removeDoublets.py \
            --input {input} \
            --output {output} \
            --prefix noDoublets
        """


rule renameSamples:
    input:
        "noDoublets.h5ad"

    output:
        "clean.h5ad"

    shell:
        """
        python src/renameSamples.py \
            --input {input} \
            --output {output}
        """


rule analysis:
    input:
        "clean.h5ad"

    output:
        "analysed.h5ad"

    shell:
        """
        python src/analysis.py \
            --input {input} \
            --output {output} \
            --prefix analysis
        """


rule diffPeaks:
    input:
        "analysed.h5ad"

    output:
        significant="{comparison}_significant_results.csv",
        all_results="{comparison}_all_results.csv"

    params:
        group1=lambda wildcards: COMPARISONS[wildcards.comparison]["group1"],
        group2=lambda wildcards: COMPARISONS[wildcards.comparison]["group2"]

    shell:
        """
        python src/diffPeaks.py \
            --input {input} \
            --output {wildcards.comparison}.h5ad \
            --group1 {params.group1} \
            --group2 {params.group2} \
            --prefix {wildcards.comparison}
        """


rule annotatePeaks:
    input:
        "{comparison}_significant_results.csv"

    output:
        "{comparison}_annotated.csv"

    shell:
        """
        python src/annotatePeaks.py \
            --peaks {input} \
            --prefix {wildcards.comparison}
        """


rule motifEnrich:
    input:
        "{comparison}_all_results.csv"

    output:
        gain="motifs_{comparison}_gain_motif_enrichment.csv",
        loss="motifs_{comparison}_loss_motif_enrichment.csv"

    shell:
        """
        python src/motifEnrich.py \
            --peaks {input} \
            --prefix motifs_{wildcards.comparison}
        """
