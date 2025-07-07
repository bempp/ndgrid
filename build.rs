#[cfg(not(feature = "scotch"))]
fn main() {
}

#[cfg(feature = "scotch")]
fn main() {
    println!("cargo:rustc-link-lib=dylib=scotch");
}
