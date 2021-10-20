//! This is a regular crate doc comment, but it also contains a partial
//! Cargo manifest.  Note the use of a *fenced* code block, and the
//! `cargo` "language".
//!
//! ```cargo
//! [dependencies]
//! rust-htslib = "0.38"
//! anyhow = "1.0"
//! ```

use anyhow::Result;
use rust_htslib::bcf::{Format, Read, Reader, Writer};

fn main() -> Result<()> {
    snakemake.redirect_stderr(&snakemake.log[0])?;
    let mut bcf_in = Reader::from_path(&snakemake.input[0])?;
    let mut updated_header = rust_htslib::bcf::header::Header::from_template(bcf_in.header());
    updated_header.push_record(b"##FORMAT=<ID=AD,Number=R,Type=Integer,Description=\"Allele depth, number of filtered reads supporting the alleles, ordered as listed in REF and ALT.\">");
    updated_header
        .push_record(b"##FORMAT=<ID=AF,Number=R,Type=Float,Description=\"Allele frequency\">");
    let mut bcf_out = Writer::from_path(&snakemake.output[0], &updated_header, true, Format::Vcf)?;
    for result in bcf_in.records() {
        let mut record = result?;
        bcf_out.translate(&mut record);
        let field_values = record.format(b"CLCAD2").integer()?;
        let allelic_fraction = field_values
            .iter()
            .map(|&allelic_depth| allelic_depth[1] as f32 / allelic_depth[0] as f32)
            .collect::<Vec<f32>>();
        record.push_format_integer(b"AD", &field_values.concat())?;
        record.push_format_float(b"AF", &allelic_fraction)?;
        bcf_out.write(&record)?;
    }
    Ok(())
}
