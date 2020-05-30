rel Edge(x : Int, y : Int)
rel Path(x : Int, y : Int)

type alias GraphPath = #{Edge, Path}

def main() : Bool =
    run(1, 0) |= Path(1, 20).

def run(amount : Int, counter : Int) : GraphPath =
    if (counter > 12 )
        #{}
    else
        let _t = (solve genEdges(amount) <+> genRules(1));
        run(amount*2, counter + 1)


def genEdges(amount : Int) : GraphPath =
    let fold = (acc : GraphPath, s:Int) -> {
            let d = s + 1;
            acc <+> Edge(s, d).
        };
        let facts = List.range(1, amount + 1) |> List.foldLeft(fold, #{});
        facts

def genRules(amount : Int) : GraphPath =
    let fold = (acc : GraphPath, _s:Int) -> {
        acc <+> #{Path(a, c) :- Path(a, b), Edge(b, c).}
    };
    let rules = List.range(1, amount + 1) |> List.foldLeft(fold, #{});
    rules <+> #{Path(a, b) :- Edge(a, b).}