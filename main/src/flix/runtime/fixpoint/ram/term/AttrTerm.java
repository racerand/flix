package flix.runtime.fixpoint.ram.term;

import flix.runtime.fixpoint.ram.RowVariable;

import java.io.PrintStream;
import java.util.Objects;

public final class AttrTerm implements RamTerm {
    private final RowVariable localVar;
    private final int index;

    public AttrTerm(RowVariable localVar, int index) {
        if (localVar == null) {
            throw new IllegalArgumentException("'localVar' must be non-null");
        }
        this.localVar = localVar;
        this.index = index;
    }

    public RowVariable getLocalVar() {
        return localVar;
    }

    public int getIndex() {
        return index;
    }

    @Override
    public void prettyPrint(PrintStream stream) {
        stream.print(toString());
    }

    @Override
    public boolean equals(Object o) {
        if (this == o) return true;
        if (o == null || getClass() != o.getClass()) return false;
        AttrTerm attrTerm = (AttrTerm) o;
        return index == attrTerm.index &&
                localVar.equals(attrTerm.localVar);
    }

    @Override
    public int hashCode() {
        return Objects.hash(localVar, index);
    }

    @Override
    public String toString() {
        return "AttrTerm{" +
                "localVar=" + localVar +
                ", index=" + index +
                '}';
    }
}