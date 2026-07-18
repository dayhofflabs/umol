// Throwaway: verify the NH3 + HCl -> NH4+ + Cl- reaction rule for a whitepaper figure.
use umol_ast::dsl::ReactionDsl;
use umol_edn::FromEdn;

const RXN: &str = r##"{:lhs {:atoms [[:N "N#h3"] [:H "H"] [:Cl "Cl"]]
        :bonds [{:id :hcl :atoms [:H :Cl] :type :single}]}
 :deltas [{:bond {:add [:N :H :single]}}
          {:bond {:remove :hcl}}
          {:atom {:modify [:N "#c+"]}}
          {:atom {:modify [:Cl "#c-"]}}]}"##;

fn main() {
    match ReactionDsl::from_edn_str(RXN) {
        Ok(rxn) => println!("PARSED OK. Canonical serialization:\n{rxn}"),
        Err(e) => {
            eprintln!("parse failed: {e}");
            std::process::exit(1);
        }
    }
}
