rel DirectedEdge(x: Int, y:Int)
rel DirectedPath(x: Int, y:Int)

def main(): Bool =
    (solve genEdges(20) <+> directedPathRule() <+> DirectedEdge(22, 23). <+> litRules() ) |= DirectedPath(1, 20).

def genEdges(x:Int): Q =
    let fold = (acc:Q, s:Int) -> {
        let d = s + 1;
        acc <+> DirectedEdge(s, d).
    };
    let facts = List.range(1, x) |> List.foldLeft(fold, #{});
    facts

def directedPathRule(): Q = #{
    DirectedPath(x, y) :- DirectedEdge(x, y).
    DirectedPath(x, z) :- DirectedPath(x, y), not DirectedEdge(y, q), DirectedEdge(q, z).
    DirectedPath(x, z) :- DirectedPath(x, y), not DirectedEdge(z, 20), DirectedEdge(y, z).
}

def litRules(): Q = #{
    DirectedPath(21, y) :- DirectedEdge(_, y).
}

type alias Q = #{DirectedEdge, DirectedPath}