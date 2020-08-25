/*
 * Copyright 2018 Magnus Madsen
 *
 * Licensed under the Apache License, Version 2.0 (the "License");
 * you may not use this file except in compliance with the License.
 * You may obtain a copy of the License at
 *
 *   http://www.apache.org/licenses/LICENSE-2.0
 *
 * Unless required by applicable law or agreed to in writing, software
 * distributed under the License is distributed on an "AS IS" BASIS,
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 * See the License for the specific language governing permissions and
 * limitations under the License.
 */

package flix.runtime.fixpoint;

import flix.runtime.fixpoint.predicate.AtomPredicate;
import flix.runtime.fixpoint.predicate.Predicate;
import flix.runtime.fixpoint.ram.interpreter.RamInterpreter;
import flix.runtime.fixpoint.ram.stmt.Stmt;
import flix.runtime.fixpoint.symbol.PredSym;

import java.util.Arrays;
import java.util.LinkedList;

public final class Solver {

    /**
     * Returns the composition of `cs1` with `cs2`.
     */
    public static ConstraintSystem compose(ConstraintSystem cs1, ConstraintSystem cs2) {
        var constraints = concat(cs1.getConstraints(), cs2.getConstraints());
        return ConstraintSystem.of(constraints);
    }

    /**
     * Returns the projection of the given predicate symbol `sym` of the given constraint system `cs`.
     */
    public static ConstraintSystem project(PredSym sym, ConstraintSystem cs) {
        if (sym == null)
            throw new IllegalArgumentException("'sym' must be non-null.");
        if (cs == null)
            throw new IllegalArgumentException("'cs' must be non-null.");

        // Collect all facts with `sym` in its head.
        var result = new LinkedList<Constraint>();
        for (Constraint fact : cs.getFacts()) {
            Predicate head = fact.getHeadPredicate();
            if (head instanceof AtomPredicate) {
                if (((AtomPredicate) head).getSym().equals(sym)) {
                    result.add(fact);
                }
            }
        }

        Constraint[] facts = result.toArray(new Constraint[0]);
        Constraint[] rules = new Constraint[0];
        Constraint[] constraints = concat(facts, rules);
        return ConstraintSystem.of(constraints);
    }

    /**
     * Solves the given constraint system `cs` with the given stratification `stf` and options `o`.
     */
    public static ConstraintSystem solve(ConstraintSystem cs, Stratification stf, Options o) {
        ca.uwaterloo.flix.runtime.solver.Solver solver = new ca.uwaterloo.flix.runtime.solver.Solver(cs, stf, o);

        long compilerMedian = evalCompiler(cs, stf, o, 0);
        long flixSolverMedian = evalFlixSolver(cs, stf, o, 0);
        long interpreterMedian = evalInterpreter(cs, stf, o, 0);
        //System.out.printf("Experiment medians, %d, %d, %d%n", compilerMedian, interpreterMedian, flixSolverMedian);

        increasingInputExperiments(cs, stf, o);
        return ConstraintSystem.of(new Constraint[0]);
    }

    private static long evalFlixSolver(ConstraintSystem cs, Stratification stf, Options o, int timesToExperiment) {
        if (timesToExperiment == 0) return 0;
        System.out.println("Start solver experiment");
        long[] times = new long[timesToExperiment];
        for (int i = 0; i < timesToExperiment; i++) {
            ca.uwaterloo.flix.runtime.solver.Solver solver = new ca.uwaterloo.flix.runtime.solver.Solver(cs, stf, o);
            long startTime = System.nanoTime();
            solver.solve();
            long endTime = System.nanoTime();
            times[i] = endTime - startTime;
            printPercentElapsed(i + 1, timesToExperiment);
        }
        long median = findMedian(times);
        long mean = findMean(times);
        System.out.printf("Median of %d runs of the Flix solver, %d%n", timesToExperiment, median);
        System.out.printf("Mean of %d runs of the Flix solver, %d%n", timesToExperiment, mean);

        return mean;
    }

    private static long evalInterpreter(ConstraintSystem cs, Stratification stf, Options o, int timesToExperiment) {
        if (timesToExperiment == 0) return 0;
        System.out.println("Start interpreter experiment");

        long[] times = new long[timesToExperiment];
        for (int i = 0; i < timesToExperiment; i++) {
            long startTime = System.nanoTime();
            Stmt compiled = DatalogCompiler.compileProgram(cs, stf, o);
            RamInterpreter.run(compiled);
            long endTime = System.nanoTime();
            times[i] = endTime - startTime;
            printPercentElapsed(i + 1, timesToExperiment);
        }
        long median = findMedian(times);
        long mean = findMean(times);
        System.out.printf("Median of %d runs of the interpreter, %d%n", timesToExperiment, median);
        System.out.printf("Mean of %d runs of the interpreter, %d%n", timesToExperiment, mean);

        return median;
    }

    private static void increasingInputExperiments(ConstraintSystem cs, Stratification stf, Options o) {
        int inputEdges = cs.getFacts().length;
        int inputRules = cs.getRules().length;

        Stmt compiled = DatalogCompiler.compileProgram(cs, stf, o);
        int programSize = compiled.size();

        System.out.printf("%d, %d, %d%n", inputEdges, inputRules, programSize);
    }

    private static long evalCompiler(ConstraintSystem cs, Stratification stf, Options o, int timesToExperiment) {
        if (timesToExperiment == 0) return 0;
        System.out.println("Start compiler vs solver experiment:");

        long[] times = new long[timesToExperiment];
        for (int i = 0; i < timesToExperiment; i++) {
            long startTime = System.nanoTime();
            Stmt compiled = DatalogCompiler.compileProgram(cs, stf, o);
            long endTime = System.nanoTime();
            times[i] = endTime - startTime;
            printPercentElapsed(i + 1, timesToExperiment);
        }
        long median = findMedian(times);
        long mean = findMean(times);
        System.out.printf("Median of %d runs of the compiler, %d%n", timesToExperiment, median);
        System.out.printf("Mean of %d runs of the compiler, %d%n", timesToExperiment, mean);


        return median;
    }

    private static void printPercentElapsed(int elapsed, int finished) {
        float percentElapsed = (float) elapsed * 100 / finished;
        if (1 >= finished / 100)
            System.out.printf("Completed %f %% of the run%n", percentElapsed);
        else if (percentElapsed % 1 < 100f / finished - 0.00001)
            System.out.printf("Completed %f %% of the run%n", percentElapsed);
    }

    private static long findMean(long[] values) {
        long sum = 0;
        for (long value : values) sum += value;
        return sum / values.length;
    }

    public static long findMedian(long[] values) {
        if (values.length == 1) return values[0];
        Arrays.sort(values);
        //System.out.println("Length = " + values.length);
        int midpoint = values.length / 2 - 1;
        //System.out.println("Midpoint = " + midpoint);
        if (values.length % 2 == 0) {
            //System.out.println("Midpoint2 = " + (midpoint + 1));
            return (values[midpoint] + values[midpoint + 1]) / 2;
        }
        return values[midpoint];
    }

    /**
     * Returns `true` if all facts in `cs2` are included in `cs1`.
     */
    public static boolean entails(ConstraintSystem cs1, ConstraintSystem cs2) {
        var entails = true;
        for (Constraint fact2 : cs2.getFacts()) {
            var found = false;
            for (Constraint fact1 : cs1.getFacts()) {
                if (fact1.entails(fact2)) {
                    found = true;
                }
            }
            if (!found) {
                entails = false;
            }
        }
        return entails;
    }

    /**
     * Returns the concatenation of the two given arrays.
     */
    private static <T> T[] concat(T[] fst, T[] snd) {
        T[] result = Arrays.copyOf(fst, fst.length + snd.length);
        System.arraycopy(snd, 0, result, fst.length, snd.length);
        return result;
    }

}
