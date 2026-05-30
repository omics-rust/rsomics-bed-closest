use std::io;
use std::path::PathBuf;

use clap::Parser;
use rsomics_common::{CommonFlags, Result, RsomicsError, Tool, ToolMeta};
use rsomics_help::{Example, FlagSpec, HelpSpec, Origin, Section};

use rsomics_bed_closest::{closest, closest_stdin};

pub const META: ToolMeta = ToolMeta {
    name: env!("CARGO_PKG_NAME"),
    version: env!("CARGO_PKG_VERSION"),
};

#[derive(Parser, Debug)]
#[command(name = "rsomics-bed-closest", disable_help_flag = true)]
pub struct Cli {
    /// First BED file A (default: stdin)
    input: Option<PathBuf>,
    /// Second BED file B (required; searched for closest features)
    #[arg(short = 'b', long)]
    b: PathBuf,
    /// Output BED (default: stdout)
    #[arg(short = 'o', long)]
    output: Option<PathBuf>,
    /// Append a signed-distance column (bedtools closest -d)
    #[arg(short = 'd', long = "distance")]
    distance: bool,
    #[command(flatten)]
    pub common: CommonFlags,
}

impl Tool for Cli {
    fn meta() -> ToolMeta {
        META
    }
    fn common(&self) -> &CommonFlags {
        &self.common
    }

    fn execute(self) -> Result<()> {
        let mut stdout_lock;
        let mut file_out;
        let out: &mut dyn io::Write = if let Some(ref p) = self.output {
            file_out = std::fs::File::create(p).map_err(RsomicsError::Io)?;
            &mut file_out
        } else {
            stdout_lock = io::stdout().lock();
            &mut stdout_lock
        };
        match self.input {
            Some(ref p) => closest(p.as_path(), &self.b, self.distance, out),
            None => closest_stdin(&self.b, self.distance, out),
        }
    }
}

pub const HELP: HelpSpec = HelpSpec {
    name: META.name,
    version: META.version,
    tagline: "Find the closest feature in B for each A interval (bedtools closest equivalent).",
    origin: Some(Origin {
        upstream: "bedtools",
        upstream_license: "MIT",
        our_license: "MIT OR Apache-2.0",
        paper_doi: Some("10.1093/bioinformatics/btq033"),
    }),
    usage_lines: &["-b <B> [OPTIONS] [INPUT]"],
    sections: &[Section {
        title: "OPTIONS",
        flags: &[
            FlagSpec {
                short: Some('b'),
                long: "b",
                aliases: &[],
                value: Some("<FILE>"),
                type_hint: Some("Path"),
                required: true,
                default: None,
                description: "Second BED file to search for closest features",
                why_default: None,
            },
            FlagSpec {
                short: Some('o'),
                long: "output",
                aliases: &[],
                value: Some("<path>"),
                type_hint: Some("Path"),
                required: false,
                default: Some("stdout"),
                description: "Output BED path",
                why_default: None,
            },
            FlagSpec {
                short: Some('d'),
                long: "distance",
                aliases: &[],
                value: None,
                type_hint: Some("bool"),
                required: false,
                default: None,
                description: "Append a signed-distance column (bedtools closest -d)",
                why_default: None,
            },
            FlagSpec {
                short: Some('h'),
                long: "help",
                aliases: &[],
                value: None,
                type_hint: Some("bool"),
                required: false,
                default: None,
                description: "Show this help",
                why_default: None,
            },
        ],
    }],
    examples: &[Example {
        description: "Find closest gene for each ATAC-seq peak",
        command: "rsomics-bed-closest peaks.bed -b genes.bed",
    }],
    json_result_schema_doc: None,
};

#[cfg(test)]
mod tests {
    use clap::CommandFactory;
    #[test]
    fn cli_definition_is_valid() {
        super::Cli::command().debug_assert();
    }
}
