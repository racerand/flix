package flix.runtime.fixpoint.ram;

import java.io.PrintStream;

public class EqualsBoolExp implements BoolExp {
    RamTerm term1;
    RamTerm term2;

    public EqualsBoolExp(RamTerm term1, RamTerm term2) {
        this.term1 = term1;
        this.term2 = term2;
    }

    @Override
    public void prettyPrint(PrintStream stream) {
        term1.prettyPrint(stream);
        stream.print(" = ");
        term2.prettyPrint(stream);
    }
}
