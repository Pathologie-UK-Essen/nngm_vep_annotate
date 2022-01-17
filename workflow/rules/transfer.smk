rule vcf2ship:
    input:
        vcf="{output_dir}/{{sample}}.annotated.filtered_ann.vcf".format(
            output_dir=config["general"]["output_path"]
        ),
        token=config["ship"]["token_path"],
    output:
        "results/transfered/{sample}.vcf",
    params:
        endpoint_url=config["ship"]["endpoint_url"],
    scripts:
        "transfer_vcf.py"