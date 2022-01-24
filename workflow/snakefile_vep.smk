import time


configfile: "config/config.yaml"


def replace_special_chars(filename):
    special_chars = [(" ", "_"), ("(", ""), (")", "")]
    for key, value in special_chars:
        filename = filename.replace(key, value)
    return filename


# TODO Limit to nNGM-Files only
samples = {
    replace_special_chars(filename[:-4]): os.path.join(root, filename)
    for root, _, files in os.walk(config["general"]["input_path"])
    for filename in files
    if (filename.endswith(".vcf") or filename.endswith(".bcf"))
    and (
        time.time() - os.path.getmtime(os.path.join(root, filename))
        > config["general"]["file_age"] * 60
    )
}


include: "rules/annotation.smk"


rule all:
    input:
        expand(
            "{output_dir}/{sample}.annotated.filtered_ann.vcf",
            output_dir=config["general"]["output_path"],
            sample=samples.keys(),
        ),
