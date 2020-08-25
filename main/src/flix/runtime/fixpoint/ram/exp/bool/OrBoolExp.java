package flix.runtime.fixpoint.ram.exp.bool;

import java.io.PrintStream;

public class OrBoolExp implements BoolExp {

    private final BoolExp leftExp;
    private final BoolExp rightExp;

    public OrBoolExp(BoolExp leftExp, BoolExp rightExp) {
        this.leftExp = leftExp;
        this.rightExp = rightExp;
    }

    @Override
    public void prettyPrint(PrintStream stream, int indentLevel) {
        leftExp.prettyPrint(stream, indentLevel);
        stream.print("||");
        rightExp.prettyPrint(stream, indentLevel);
    }

    public BoolExp getLeftExp() {
        return leftExp;
    }

    public BoolExp getRightExp() {
        return rightExp;
    }

    @Override
    public String toString() {
        return "OrBoolExp{" +
                "leftExp=" + leftExp +
                ", rightExp=" + rightExp +
                '}';
    }
}