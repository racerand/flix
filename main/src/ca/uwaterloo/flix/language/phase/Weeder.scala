package ca.uwaterloo.flix.language.phase

import ca.uwaterloo.flix.language.Compiler
import ca.uwaterloo.flix.language.ast.{WeededAst, SourceLocation, ParsedAst, Name}
import ca.uwaterloo.flix.util.Validation
import Validation._

import scala.collection.mutable

/**
 * The Weeder phase performs simple syntactic checks and rewritings.
 */
// TODO: JoinSemiLattice vs. CompleteLattice.
// TODO: valid "Traits"

object Weeder {

  import WeederError._

  /**
   * A common super-type for weeding errors.
   */
  sealed trait WeederError extends Compiler.CompilationError

  object WeederError {

    implicit val consoleCtx = Compiler.ConsoleCtx

    /**
     * An error raised to indicate that the attribute `name` was declared multiple times.
     *
     * @param name the name of the attribute.
     * @param location1 the location of the first declaration.
     * @param location2 the location of the second declaration.
     */
    case class DuplicateAttribute(name: String, location1: SourceLocation, location2: SourceLocation) extends WeederError {
      val format =
        s"""${consoleCtx.blue(s"-- SYNTAX ERROR ------------------------------------------------ ${location1.formatSource}")}
            |
            |${consoleCtx.red(s">> Duplicate attribute name '$name'.")}
            |
            |First declaration was here:
            |${location1.underline}
            |Second declaration was here:
            |${location2.underline}
         """.stripMargin
    }

    /**
     * An error raised to indicate that the formal argument `name` was declared multiple times.
     *
     * @param name the name of the argument.
     * @param location1 the location of the first declaration.
     * @param location2 the location of the second declaration.
     */
    case class DuplicateFormal(name: String, location1: SourceLocation, location2: SourceLocation) extends WeederError {
      val format =
        s"""${consoleCtx.blue(s"-- SYNTAX ERROR ------------------------------------------------ ${location1.formatSource}")}
            |
            |${consoleCtx.red(s">> Duplicate formal argument '$name'.")}
            |
            |First declaration was here:
            |${location1.underline}
            |Second declaration was here:
            |${location2.underline}
         """.stripMargin
    }

    /**
     * An error raised to indicate that the tag `name` was declared multiple times.
     *
     * @param name the name of the tag.
     * @param location1 the location of the first declaration.
     * @param location2 the location of the second declaration.
     */
    case class DuplicateTag(name: String, location1: SourceLocation, location2: SourceLocation) extends WeederError {
      val format =
        s"""${consoleCtx.blue(s"-- SYNTAX ERROR ------------------------------------------------ ${location1.formatSource}")}
            |
            |${consoleCtx.red(s">> Duplicate tag name '$name'.")}
            |
            |First declaration was here:
            |${location1.underline}
            |Second declaration was here:
            |${location2.underline}
         """.stripMargin
    }

    /**
     * An error raised to indicate that an illegal term occurs inside a predicate.
     *
     * @param location the location where the illegal term occurs.
     */
    case class IllegalTerm(message: String, location: SourceLocation) extends WeederError {
      val format = s"Error: $message Term was here: ${location.format}.\n"
    }

    /**
     * An error raised to indicate an illegal lattice definition.
     *
     * @param location the location where the illegal definition occurs.
     */
    case class IllegalLattice(message: String, location: SourceLocation) extends WeederError {
      val format = s"Error: $message Definition was here: ${location.format}.\n"
    }

    /**
     * An error raised to indicate that the variable `name` occurs multiple times in the same pattern.
     *
     * @param name the name of the variable.
     *
     * @param location1 the location of the first use of the variable.
     * @param location2 the location of the second use of the variable.
     */
    case class NonLinearPattern(name: String, location1: SourceLocation, location2: SourceLocation) extends WeederError {
      val format =
        s"Error: Non-linear pattern: The variable: '$name' occurs twice.\n" +
          s"  First occurrence was here: ${location1.format}\n" +
          s"  Second occurrence was here: ${location2.format}\n"
    }

    /**
     * An error raised to indicate that a syntactic construct, although successfully parsed, is currently not supported.
     *
     * @param message the error message.
     * @param location the location of the syntactic construct.
     */
    case class Unsupported(message: String, location: SourceLocation) extends WeederError {
      val format =
        s"Error: $message at ${location.format}\n."
    }

  }

  /**
   * Compiles the given parsed `past` to a weeded ast.
   */
  def weed(past: ParsedAst.Root): Validation[WeededAst.Root, WeederError] = {
    @@(past.declarations.map(Declaration.compile)) map WeededAst.Root
  }

  object Declaration {

    /**
     * Compiles the given parsed declaration `past` to a weeded declaration.
     */
    def compile(past: ParsedAst.Declaration): Validation[WeededAst.Declaration, WeederError] = past match {
      case d: ParsedAst.Declaration.Namespace => compile(d)
      case d: ParsedAst.Declaration.Fact => compile(d)
      case d: ParsedAst.Declaration.Rule => compile(d)
      case d: ParsedAst.Definition => Definition.compile(d)
    }

    /**
     * Compiles the given parsed namespace declaration `past` to a weeded namespace declaration.
     */
    def compile(past: ParsedAst.Declaration.Namespace): Validation[WeededAst.Declaration.Namespace, WeederError] =
      @@(past.body.map(compile)) map (ds => WeededAst.Declaration.Namespace(past.name, ds))

    /**
     * Compiles the given parsed fact `past` to a weeded fact.
     */
    def compile(past: ParsedAst.Declaration.Fact): Validation[WeededAst.Declaration.Fact, WeederError] =
      Predicate.Head.compile(past.head) map {
        case p => WeededAst.Declaration.Fact(p)
      }

    /**
     * Compiles the parsed rule `past` to a weeded rule.
     */
    def compile(past: ParsedAst.Declaration.Rule): Validation[WeededAst.Declaration.Rule, WeederError] = {
      val headVal = Predicate.Head.compile(past.head)
      val bodyVal = @@(past.body.map(Predicate.Body.compile))

      @@(headVal, bodyVal) map {
        case (head, body) => WeededAst.Declaration.Rule(head, body)
      }
    }

  }

  object Definition {

    /**
     * Compiles the given parsed definition `past` to a weeded definition.
     */
    def compile(past: ParsedAst.Definition): Validation[WeededAst.Declaration, WeederError] = past match {
      case d: ParsedAst.Definition.Value => compile(d)
      case d: ParsedAst.Definition.Function => compile(d)
      case d: ParsedAst.Definition.Enum => compile(d)
      case d: ParsedAst.Definition.Lattice => compile(d)
      case d: ParsedAst.Definition.Relation => compile(d)
    }

    /**
     * Compiles the given parsed value declaration `past` to a weeded definition.
     */
    def compile(past: ParsedAst.Definition.Value): Validation[WeededAst.Definition.Constant, WeederError] =
      @@(Expression.compile(past.e), Type.compile(past.tpe)) map {
        case (exp, tpe) => WeededAst.Definition.Constant(past.ident, exp, tpe, past.loc)
      }

    /**
     * Compiles the given parsed function declaration `past` to a weeded definition.
     */
    def compile(past: ParsedAst.Definition.Function): Validation[WeededAst.Definition.Constant, WeederError] = {
      val seen = mutable.Map.empty[String, Name.Ident]

      val formalsVal = @@(past.formals.map {
        case ParsedAst.FormalArg(ident, tpe) => seen.get(ident.name) match {
          case None =>
            seen += (ident.name -> ident)
            Type.compile(tpe) map (t => WeededAst.FormalArg(ident, t))
          case Some(otherIdent) =>
            (DuplicateFormal(ident.name, otherIdent.location, ident.location): WeederError).toFailure
        }
      })

      @@(formalsVal, Expression.compile(past.body), Type.compile(past.tpe)) map {
        case (args, body, retType) =>
          val exp = WeededAst.Expression.Lambda(args, body, retType, past.body.loc)
          val tpe = WeededAst.Type.Function(args map (_.tpe), retType)
          WeededAst.Definition.Constant(past.ident, exp, tpe, past.loc)
      }
    }

    /**
     * Compiles the given parsed enum declaration `past` to a weeded enum definition.
     *
     * Returns [[Failure]] if the same tag name occurs twice.
     */
    def compile(past: ParsedAst.Definition.Enum): Validation[WeededAst.Definition.Enum, WeederError] =
      Validation.fold[ParsedAst.Type.Tag, Map[String, WeededAst.Type.Tag], WeederError](past.cases, Map.empty) {
        // loop through each tag declaration.
        case (macc, tag@ParsedAst.Type.Tag(Name.Ident(name, location), _)) => macc.get(name) match {
          // check if the tag was already declared.
          case None => Type.compile(tag) map (tpe => macc + (name -> tpe.asInstanceOf[WeededAst.Type.Tag]))
          case Some(otherTag) => DuplicateTag(name, otherTag.tagName.location, location).toFailure
        }
      } map {
        case m => WeededAst.Definition.Enum(past.ident, m, past.loc)
      }

    /**
     * Compiles the given parsed lattice `past` to a weeded lattice definition.
     */
    def compile(past: ParsedAst.Definition.Lattice): Validation[WeededAst.Definition.Lattice, WeederError] = {
      val elmsVal = @@(past.elms.toList.map(Expression.compile))
      elmsVal flatMap {
        case bot :: leq :: lub :: Nil => WeededAst.Definition.Lattice(past.ident, bot, leq, lub, past.loc).toSuccess
        case _ => IllegalLattice("Lattice definition must have exactly three elements: bot, leq and lub.", past.loc).toFailure
      }
    }

    /**
     * Compiles the given parsed relation `past` to a weeded relation definition.
     */
    def compile(past: ParsedAst.Definition.Relation): Validation[WeededAst.Definition.Relation, WeederError] = {
      val seen = mutable.Map.empty[String, Name.Ident]

      val attributesVal = past.attributes.map {
        case ParsedAst.Attribute(ident, ptype) => seen.get(ident.name) match {
          // check if the attribute name was already declared.
          case None =>
            seen += (ident.name -> ident)
            Type.compile(ptype) map (tpe => WeededAst.Attribute(ident, tpe))
          case Some(otherIdent) =>
            (DuplicateAttribute(ident.name, otherIdent.location, ident.location): WeederError).toFailure
        }
      }

      @@(attributesVal) map {
        case attributes => WeededAst.Definition.Relation(past.ident, attributes, past.loc)
      }
    }
  }

  object Literal {
    /**
     * Compiles the parsed literal `past` to a weeded literal.
     */
    def compile(past: ParsedAst.Literal): Validation[WeededAst.Literal, WeederError] = past match {
      case ParsedAst.Literal.Unit(loc) => WeededAst.Literal.Unit(loc).toSuccess
      case ParsedAst.Literal.Bool(loc, b) => WeededAst.Literal.Bool(b, loc).toSuccess
      case ParsedAst.Literal.Int(loc, i) => WeededAst.Literal.Int(i, loc).toSuccess
      case ParsedAst.Literal.Str(loc, s) => WeededAst.Literal.Str(s, loc).toSuccess
      case ParsedAst.Literal.Tag(loc, name, ident, literal) => compile(literal) map (l => WeededAst.Literal.Tag(name, ident, l, loc))
      case ParsedAst.Literal.Tuple(loc, pelms) => @@(pelms map compile) map {
        case elms => WeededAst.Literal.Tuple(elms, loc)
      }
    }
  }

  object Expression {
    /**
     * Compiles the parsed expression `past` to a weeded expression.
     */
    def compile(past: ParsedAst.Expression): Validation[WeededAst.Expression, WeederError] = past match {
      case ParsedAst.Expression.Lit(loc, literal) =>
        Literal.compile(literal) map {
          case lit => WeededAst.Expression.Lit(lit, loc)
        }

      case ParsedAst.Expression.Var(loc, name) =>
        WeededAst.Expression.Var(name, loc).toSuccess

      case ParsedAst.Expression.Apply(loc, wlambda, wargs) =>
        @@(compile(wlambda), @@(wargs map compile)) map {
          case (lambda, args) => WeededAst.Expression.Apply(lambda, args, loc)
        }

      case ParsedAst.Expression.Lambda(loc, pargs, ptype, pbody) =>
        val argsVal = @@(pargs map {
          case ParsedAst.FormalArg(ident, tpe) => Type.compile(tpe) map (t => WeededAst.FormalArg(ident, t))
        })
        @@(argsVal, Type.compile(ptype), compile(pbody)) map {
          case (args, tpe, body) => WeededAst.Expression.Lambda(args, body, tpe, loc)
        }

      case ParsedAst.Expression.Unary(loc, op, pe) =>
        compile(pe) map {
          case e => WeededAst.Expression.Unary(op, e, loc)
        }

      case ParsedAst.Expression.Binary(pe1, loc, op, pe2) =>
        @@(compile(pe1), compile(pe2)) map {
          case (e1, e2) => WeededAst.Expression.Binary(op, e1, e2, loc)
        }

      case ParsedAst.Expression.IfThenElse(loc, pe1, pe2, pe3) =>
        @@(compile(pe1), compile(pe2), compile(pe3)) map {
          case (e1, e2, e3) => WeededAst.Expression.IfThenElse(e1, e2, e3, loc)
        }

      case ParsedAst.Expression.Let(loc, ident, pvalue, pbody) =>
        @@(compile(pvalue), compile(pbody)) map {
          case (value, body) => WeededAst.Expression.Let(ident, value, body, loc)
        }

      case ParsedAst.Expression.Match(loc, pe, prules) =>
        val rulesVal = prules map {
          case (pat, body) => @@(Pattern.compile(pat), compile(body))
        }
        @@(compile(pe), @@(rulesVal)) map {
          case (e, rs) => WeededAst.Expression.Match(e, rs, loc)
        }

      case ParsedAst.Expression.Infix(pe1, loc, name, pe2) =>
        @@(compile(pe1), compile(pe2)) map {
          case (e1, e2) => WeededAst.Expression.Apply(WeededAst.Expression.Var(name, loc), List(e1, e2), loc)
        }

      case ParsedAst.Expression.Tag(loc, enumName, tagName, pe) => compile(pe) map {
        case e => WeededAst.Expression.Tag(enumName, tagName, e, loc)
      }

      case ParsedAst.Expression.Tuple(loc, pelms) => @@(pelms map compile) map {
        case elms => WeededAst.Expression.Tuple(elms, loc)
      }

      case ParsedAst.Expression.Ascribe(loc, pe, ptype) =>
        @@(compile(pe), Type.compile(ptype)) map {
          case (e, tpe) => WeededAst.Expression.Ascribe(e, tpe, loc)
        }

      case ParsedAst.Expression.Error(loc, ptype) =>
        Type.compile(ptype) map {
          case tpe => WeededAst.Expression.Error(tpe, loc)
        }
    }
  }

  object Pattern {
    /**
     * Compiles the parsed pattern `past`.
     *
     * Returns [[Failure]] if the pattern is non-linear, i.e. if the same variable occurs twice.
     */
    def compile(past: ParsedAst.Pattern): Validation[WeededAst.Pattern, WeederError] = {
      val seen = mutable.Map.empty[String, Name.Ident]

      def visit(p: ParsedAst.Pattern): Validation[WeededAst.Pattern, WeederError] = p match {
        case ParsedAst.Pattern.Wildcard(loc) => WeededAst.Pattern.Wildcard(loc).toSuccess
        case ParsedAst.Pattern.Var(loc, ident) => seen.get(ident.name) match {
          case None =>
            seen += (ident.name -> ident)
            WeededAst.Pattern.Var(ident, loc).toSuccess
          case Some(otherIdent) =>
            NonLinearPattern(ident.name, otherIdent.location, ident.location).toFailure
        }
        case ParsedAst.Pattern.Lit(loc, literal) => Literal.compile(literal) map {
          case lit => WeededAst.Pattern.Lit(lit, loc)
        }
        case ParsedAst.Pattern.Tag(loc, enumName, tagName, ppat) => visit(ppat) map {
          case pat => WeededAst.Pattern.Tag(enumName, tagName, pat, loc)
        }
        case ParsedAst.Pattern.Tuple(loc, pelms) => @@(pelms map visit) map {
          case elms => WeededAst.Pattern.Tuple(elms, loc)
        }
      }

      visit(past)
    }
  }


  object Predicate {

    object Head {

      /**
       * Compiles the given parsed predicate `p` to a weeded predicate.
       */
      def compile(past: ParsedAst.Predicate): Validation[WeededAst.Predicate.Head, WeederError] = past match {
        case p: ParsedAst.Predicate.Unresolved => compile(p)
      }

      /**
       * Compiles the given parsed predicate `past` to a weeded predicate.
       */
      def compile(past: ParsedAst.Predicate.Unresolved): Validation[WeededAst.Predicate.Head, WeederError] =
        @@(past.terms.map(Term.Head.compile)) map {
          case terms => WeededAst.Predicate.Head(past.name, terms, past.loc)
        }
    }

    object Body {

      /**
       * Compiles the given parsed predicate `p` to a weeded predicate.
       */
      def compile(past: ParsedAst.Predicate): Validation[WeededAst.Predicate.Body, WeederError] = past match {
        case p: ParsedAst.Predicate.Unresolved => compile(p)
      }

      /**
       * Compiles the given parsed predicate `p` to a weeded predicate.
       */
      def compile(past: ParsedAst.Predicate.Unresolved): Validation[WeededAst.Predicate.Body, WeederError] =
        @@(past.terms.map(Term.Body.compile)) map {
          case terms => WeededAst.Predicate.Body(past.name, terms, past.loc)
        }
    }

  }

  object Term {

    object Head {

      /**
       * Compiles the given parsed head term `past` to a weeded term.
       *
       * Returns [[Failure]] if the term contains a wildcard variable.
       */
      def compile(past: ParsedAst.Term): Validation[WeededAst.Term.Head, WeederError] = past match {
        case ParsedAst.Term.Wildcard(loc) => IllegalTerm("Wildcard variables not allowed in head terms.", loc).toFailure
        case ParsedAst.Term.Var(loc, ident) => WeededAst.Term.Head.Var(ident, loc).toSuccess
        case ParsedAst.Term.Lit(loc, literal) => Literal.compile(literal) map {
          case lit => WeededAst.Term.Head.Lit(lit, loc)
        }
        case ParsedAst.Term.Ascribe(loc, pterm, ptpe) =>
          @@(compile(pterm), Type.compile(ptpe)) map {
            case (term, tpe) => WeededAst.Term.Head.Ascribe(term, tpe, loc)
          }
        case ParsedAst.Term.Apply(loc, name, pargs) =>
          @@(pargs map compile) map {
            case args => WeededAst.Term.Head.Apply(name, args, loc)
          }
      }
    }

    object Body {
      /**
       * Compiles the given parsed body term `past` to a weeded term.
       *
       * Returns [[Failure]] if the term contains a function call.
       */
      def compile(past: ParsedAst.Term): Validation[WeededAst.Term.Body, WeederError] = past match {
        case ParsedAst.Term.Wildcard(loc) => WeededAst.Term.Body.Wildcard(loc).toSuccess
        case ParsedAst.Term.Var(loc, ident) => WeededAst.Term.Body.Var(ident, loc).toSuccess
        case ParsedAst.Term.Lit(loc, literal) => Literal.compile(literal) map {
          case lit => WeededAst.Term.Body.Lit(lit, loc)
        }
        case ParsedAst.Term.Ascribe(loc, pterm, ptpe) =>
          @@(compile(pterm), Type.compile(ptpe)) map {
            case (term, tpe) => WeededAst.Term.Body.Ascribe(term, tpe, loc)
          }
        case ParsedAst.Term.Apply(loc, name, args) => IllegalTerm("Function calls not allowed in body terms.", loc).toFailure
      }
    }

  }

  object Type {
    /**
     * Weeds the given parsed type `past`.
     */
    def compile(past: ParsedAst.Type): Validation[WeededAst.Type, WeederError] = past match {
      case ParsedAst.Type.Unit => WeededAst.Type.Unit.toSuccess
      case ParsedAst.Type.Var(name) => WeededAst.Type.Ambiguous(name).toSuccess
      case ParsedAst.Type.Function(pformals, pret) =>
        val formalTypeVal = @@(pformals map compile)
        val returnTypeVal = compile(pret)
        @@(formalTypeVal, returnTypeVal) map {
          case (formals, retTpe) => WeededAst.Type.Function(formals, retTpe)
        }
      case ParsedAst.Type.Tag(ident, ptype) => compile(ptype) map {
        case tpe => WeededAst.Type.Tag(ident, tpe)
      }
      case ParsedAst.Type.Tuple(pelms) => @@(pelms map compile) map {
        case elms => WeededAst.Type.Tuple(elms)
      }
      case ParsedAst.Type.Parametric(name, pelms) =>
        Unsupported("Parametric types are not yet supported.", name.location).toFailure
    }
  }


}
