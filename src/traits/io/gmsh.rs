//! Gmsh I/O
use crate::traits::{Builder, Grid};
use std::fs::{self, File};
use std::io::{BufRead, BufReader, Read};

pub trait GmshExport: Grid {
    //! Grid export for Gmsh

    /// Generate the Gmsh string for a grid
    fn to_gmsh_string(&self) -> String;

    /// Export as Gmsh
    fn export_as_gmsh(&self, filename: &str) {
        let gmsh_s = self.to_gmsh_string();
        fs::write(filename, gmsh_s).expect("Unable to write file");
    }
}

pub trait GmshImport: Builder {
    //! Grid import for Gmsh

    /// Generate grid from a Gmsh v1 string
    fn import_from_gmsh_v1(&mut self, s: String);

    /// Generate grid from a Gmsh v2 string
    fn import_from_gmsh_string_v2(&mut self, s: String);

    /// Generate grid from a Gmsh v2 binary
    fn import_from_gmsh_binary_v2(
        &mut self,
        reader: BufReader<File>,
        data_size: usize,
        // check endianness.
        is_le: bool,
    );

    /// Generate grid from a Gmsh v4 string
    fn import_from_gmsh_string_v4(&mut self, s: String);

    /// Generate grid from a Gmsh v4 binary
    fn import_from_gmsh_binary_v4(
        &mut self,
        reader: BufReader<File>,
        data_size: usize,
        // check endianness.
        is_le: bool,
    );

    /// Generate grid from Gmsh
    fn import_from_gmsh(&mut self, filename: &str) {
        let f = File::open(filename).expect("Unable to open file");
        let mut reader = BufReader::new(f);

        let mut line = String::new();
        reader.read_line(&mut line).expect("Unable to read header");

        if line.starts_with("$NOD") {
            let mut content = String::new();
            reader
                .read_to_string(&mut content)
                .expect("Unable to read mesh");

            content.replace_range(0..0, "$Nodes\n");
            for (f, r) in [
                ("$ENDNOD", "$EndNodes"),
                ("$ELM", "$Elements"),
                ("$ENDELM", "$EndElements"),
            ]
            .iter()
            {
                let Some(offset) = content.find(f) else {
                    panic!("Invalid file format");
                };
                content.replace_range(offset..offset + f.len(), r);
            }

            self.import_from_gmsh_v1(content);
            return;
        }

        reader.read_line(&mut line).expect("Unable to read header");

        let [version, binary_mode, data_size] = line.split("\n").collect::<Vec<_>>()[1]
            .split(" ")
            .collect::<Vec<_>>()[..]
        else {
            panic!("Unrecognised format");
        };

        const GMSH_INT_SIZE: usize = 4;

        if binary_mode == "1" {
            let data_size = data_size
                .parse::<usize>()
                .expect("Unable to parse data size");

            assert!(
                data_size <= std::mem::size_of::<usize>(),
                "Unsupported data size"
            );

            let mut buf = [0u8; GMSH_INT_SIZE];
            reader
                .read_exact(&mut buf)
                .expect("Unable to read endianness");

            let is_le = u32::from_le_bytes(buf) == 1;

            if version.starts_with("2") {
                self.import_from_gmsh_binary_v2(reader, data_size, is_le);
            } else {
                self.import_from_gmsh_binary_v4(reader, data_size, is_le);
            }
        } else {
            let mut content = String::new();
            reader
                .read_to_string(&mut content)
                .expect("Unable to read content");

            if version.starts_with("2") {
                self.import_from_gmsh_string_v2(content);
            } else {
                self.import_from_gmsh_string_v4(content);
            }
        }
    }
}
