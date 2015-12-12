package ca.uwaterloo.flix.language.backend.phase

import ca.uwaterloo.flix.language.ast.TypedAst.Expression.Var
import ca.uwaterloo.flix.language.ast.{SourcePosition, Name, SourceLocation, TypedAst}

object Verifier {

  abstract class Property(tpe: TypedAst.Type) {
    val x = mkVar("x", tpe)
    val y = mkVar("y", tpe)
    val z = mkVar("y", tpe)

    val property: Formula
  }

  object Theorems {

    case class Commutativity(op: TypedAst.Expression) {
      val (x, y) = (mkVar("x", op.tpe), mkVar("y", op.tpe))
      val f: Formula.Lambda = ???

      val property = ∀(x, y)(f(x, y) ≡ f(y, x))
    }

    /**
      * Properties of Partial Orders.
      */
    object PartialOrder {

      /**
        * Reflexivity.
        */
      case class Reflexivity(lattice: TypedAst.Definition.BoundedLattice) extends Property(lattice.tpe) {
        val property = ∀(x)(x ⊑ x)
      }

      /**
        * Anti-symmetry.
        */
      case class AntiSymmetry(lattice: TypedAst.Definition.BoundedLattice) extends Property(lattice.tpe) {
        val property = ∀(x, y)((x ⊑ y) ∧ (y ⊑ x)) → (x ≡ y)
      }

      /**
        * Transitivity.
        */
      case class Transitivity(lattice: TypedAst.Definition.BoundedLattice) extends Property(lattice.tpe) {
        val property = ∀(x, y, z)(((x ⊑ y) ∧ (y ⊑ z)) → (x ⊑ z))
      }

    }

    /**
      * Properties of Join Semi Lattices.
      */
    object JoinSemiLattice {

    }

  }


  /**
    * Returns a typed AST node that represents a variable of the given `name` and with the given type `tpe`.
    */
  private def mkVar(name: String, tpe: TypedAst.Type): Formula.Var =
    ???

  sealed trait Formula

  object Formula {

    // TODO: probably need to carry a type.

    case class BinOp(name: String) extends Formula {
      def apply(f1: Formula, f2: Formula): Formula = ???

    }

    case class Var(name: String) extends Formula

    case class Lambda() extends Formula {
      def apply(that: Formula): Formula = ???

      def apply(e1: Formula, e2: Formula): Formula = ???
    }

    case class Implication(f1: Formula, f2: Formula) extends Formula

    case class Leq(f1: Formula, f2: Formula) extends Formula

    case class Lub(f1: Formula, f2: Formula) extends Formula

  }

  def ∀(x: Formula.Var*)(f: Formula): Formula = ???

  implicit class RichFormula(val thiz: Formula) {


    def ∧(that: Formula): Formula = ???

    def →(that: Formula): Formula = Formula.Implication(thiz, that)

    def ⊑(that: Formula): Formula = Formula.Leq(thiz, that)

    def ⊔(that: Formula): Formula = Formula.Lub(thiz, that)

    def ≡(that: Formula): Formula = ???

  }

  //  /**
  //   * Least Element: ?x. ? ? x.
  //   */
  //  def leastElement(leq: Term.Abs, bot: Value): Term.Abs =
  //    Term.Abs('x, leq.typ,
  //      leq.call(bot.toTerm, 'x))
  //
  //
  //  /**
  //   * Upper Bound: ?x, y. x ? (x ? y) ? y ? (x ? y).
  //   */
  //  def upperBound(leq: Term.Abs, lub: Term.Abs): Term.Abs =
  //    Term.Abs('x, leq.typ, Term.Abs('y, leq.typ,
  //      leq.call('x, lub.call('x, 'y)) && leq.call('y, lub.call('x, 'y))))
  //
  //  /**
  //   * Least Upper Bound: ?x, y, z. x ? z ? y ? z ? x ? y ? z.
  //   */
  //  def leastUpperBound(leq: Term.Abs, lub: Term.Abs): Term.Abs =
  //    Term.Abs('x, leq.typ, Term.Abs('y, leq.typ, Term.Abs('z, leq.typ,
  //      (leq.call('x, 'z) && leq.call('y, 'z)) ==> leq.call(lub.call('x, 'y), 'z))))
  //
  //
  //  /**
  //   * Non-Negative: ?x. f(x) > 0.
  //   */
  //  def nonNegative(h: Term.Abs): Term.Abs =
  //    Term.Abs('x, h.typ,
  //      Term.BinaryOp(BinaryOperator.Greater, h.call('x), Term.Int(0)))
  //
  //  /**
  //   * Stricly-Decreasing: ?x, y. x ? y ? x != y ? f(x) > f(y).
  //   */
  //  def strictlyDecreasing(h: Term.Abs, leq: Term.Abs): Term.Abs =
  //    Term.Abs('x, h.typ, Term.Abs('y, h.typ,
  //      (leq.call('x, 'y) && (Term.Var('x) !== Term.Var('y))) ==>
  //        Term.BinaryOp(BinaryOperator.Greater, h.call('x), h.call('y))))
  //
  //  /**
  //   * Monotone: ?x1, x2. x1 ? x2 ? f(x1) ? f(x2).
  //   */
  //  def monotone1(f: Term.Abs, leq: Term.Abs): Term.Abs =
  //    Term.Abs('x1, leq.typ, Term.Abs('x2, leq.typ,
  //      leq.call('x1, 'x2) ==> leq.call(f.call('x1), f.call('x2))))
  //
  //  /**
  //   * Monotone: ?x1, x2, y1, y2. x1 ? x2 ? y1 ? y2 ? f(x1, y1) ? f(x2, y2).
  //   */
  //  def monotone2(f: Term.Abs, leq: Term.Abs): Term.Abs =
  //    Term.Abs('x1, leq.typ, Term.Abs('x2, leq.typ, Term.Abs('y1, leq.typ, Term.Abs('y2, leq.typ,
  //      (leq.call('x1, 'x2) && leq.call('y1, 'y2)) ==> leq.call(f.call('x1, 'y1), f.call('x2, 'y2))))))

}