use trybuild::TestCases;

#[test]
fn ui() {
    TestCases::new().compile_fail("tests/ui/*.rs");
}
