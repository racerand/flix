rel Edge(x : Int, y : Int)
rel Path(x : Int, y : Int)

type alias GraphPath = #{Edge, Path}

def main() : Bool =
    run(1, 0) |= Path(1, 20).

def run(amount : Int, counter : Int) : GraphPath =
    if (counter > 12 )
        #{}
    else
        let t = (solve genEdges(0) <+> genRules(amount));
        let _s = run(amount*2, counter + 1);
        t


def genEdges(amount : Int) : GraphPath =
    let fold = (acc : GraphPath, s:Int) -> {
            let d = s + 1;
            acc <+> Edge(s, d).
        };
        let facts = List.range(1, amount + 1) |> List.foldLeft(fold, #{});
        facts

def genRules(amount : Int) : GraphPath =
    let fold = (acc : GraphPath, _s:Int) -> {
        acc <+> #{Path(a, k) :- Path(a, k)/*, Edge(b, c), Edge(c, d),
            Edge(d, e), Edge(e, f), Edge(f, g), Edge(g, h), Edge(h, i), Edge(i, j), Edge(j, k)*/.}
    };
    let rules = List.range(1, amount + 1) |> List.foldLeft(fold, #{});
    rules// <+> #{Path(a, b) :- Edge(a, b).}