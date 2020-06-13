rel DirectedEdge(x: Int, y:Int)
rel DirectedPath(x: Int, y:Int)
rel NotPath(x: Int, y:Int)

type alias Q = #{DirectedEdge, DirectedPath, NotPath}

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

def notPathRules() : Q = #{
    NotPath(y, x) :- DirectedPath(x, y), not DirectedPath(y, x).
}

def main(): #{DirectedPath} =
    project DirectedPath (solve genEdges(5) <+> directedPathRule() <+> notPathRules())