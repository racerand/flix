rel DirectedEdge(x: Int, y:Int)
rel DirectedPath(x: Int, y:Int)

type alias Q = #{DirectedEdge, DirectedPath}

def genEdges(x:Int): Q =
    let fold = (acc:Q, s:Int) -> {
        let d = s + 1;
        acc <+> DirectedEdge(s, d).
    };
    let facts = List.range(1, x) |> List.foldLeft(fold, #{});
    facts

def directedPathRule(): Q = #{
    DirectedPath(x, y) :- DirectedEdge(x, y).
    DirectedPath(x, z) :- DirectedPath(x, y), DirectedEdge(y, z).
}

def main():Bool =
    (solve genEdges(20) <+> directedPathRule()) |= DirectedPath(1,20).