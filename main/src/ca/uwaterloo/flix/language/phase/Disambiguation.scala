/*
 *  Copyright 2016 Magnus Madsen
 *
 *  Licensed under the Apache License, Version 2.0 (the "License");
 *  you may not use this file except in compliance with the License.
 *  You may obtain a copy of the License at
 *
 *  http://www.apache.org/licenses/LICENSE-2.0
 *
 *  Unless required by applicable law or agreed to in writing, software
 *  distributed under the License is distributed on an "AS IS" BASIS,
 *  WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 *  See the License for the specific language governing permissions and
 *  limitations under the License.
 */

package ca.uwaterloo.flix.language.phase

import ca.uwaterloo.flix.language.ast.NamedAst.Program
import ca.uwaterloo.flix.language.ast.{Ast, Name, NamedAst}
import ca.uwaterloo.flix.language.errors.TypeError
import ca.uwaterloo.flix.language.errors.TypeError.UnresolvedDefinition
import ca.uwaterloo.flix.language.phase.Unification._
import ca.uwaterloo.flix.util.Result
import ca.uwaterloo.flix.util.Result.{Err, Ok}

import scala.collection.mutable

object Disambiguation {

  /**
    * An ADT for the result of looking up a name.
    */
  sealed trait LookupResult

  object LookupResult {

    case class Defn(ns: Name.NName, defn: NamedAst.Declaration.Definition) extends LookupResult

    case class Hook(hook: Ast.Hook) extends LookupResult

  }

  /**
    * Finds the definition with the qualified name `qname` in the namespace `ns0`.
    */
  // TODO: Better to use Result[]
  def lookupRef(qname: Name.QName, ns0: Name.NName, program: Program): InferMonad[LookupResult] = {
    // check whether the reference is fully-qualified.
    if (qname.isUnqualified) {
      // Case 1: Unqualified reference. Lookup both the definition and the hook.
      val defnOpt = program.definitions.getOrElse(ns0, Map.empty).get(qname.ident.name)
      val hookOpt = program.hooks.getOrElse(ns0, Map.empty).get(qname.ident.name)

      (defnOpt, hookOpt) match {
        case (None, None) => failM(UnresolvedDefinition(qname, ns0, qname.loc))
        case (Some(defn), None) => liftM(LookupResult.Defn(ns0, defn))
        case (None, Some(hook)) => liftM(LookupResult.Hook(hook))
        case (Some(defn), Some(hook)) => ??? // TODO: Overloaded name.
      }
    } else {
      // Case 2: Qualified. Lookup both the definition and the hook.
      val defnOpt = program.definitions.getOrElse(qname.namespace, Map.empty).get(qname.ident.name)
      val hookOpt = program.hooks.getOrElse(qname.namespace, Map.empty).get(qname.ident.name)

      (defnOpt, hookOpt) match {
        case (None, None) => failM(UnresolvedDefinition(qname, ns0, qname.loc))
        case (Some(defn), None) => liftM(LookupResult.Defn(qname.namespace, defn))
        case (None, Some(hook)) => liftM(LookupResult.Hook(hook))
        case (Some(defn), Some(hook)) => ??? // TODO: Overloaded name.
      }
    }
  }

  /**
    * Finds the enum definition matching the given qualified name and tag.
    */
  def lookupEnumByTag(qname: Name.QName, tag: Name.Ident, ns: Name.NName, program: Program): Result[NamedAst.Declaration.Enum, TypeError] = {
    /*
     * Lookup the tag name in all enums across all namespaces.
     */
    val globalMatches = mutable.Set.empty[NamedAst.Declaration.Enum]
    for ((_, decls) <- program.enums) {
      for ((enumName, decl) <- decls) {
        for ((tagName, caze) <- decl.cases) {
          if (tag.name == tagName) {
            globalMatches += decl
          }
        }
      }
    }


    // Case 1: Exact match found. Simply return it.
    if (globalMatches.size == 1) {
      return Ok(globalMatches.head)
    }

    // Case 2: No or multiple matches found.
    // Lookup the tag in either the fully qualified namespace or the current namespace.
    val namespace = if (qname.isQualified) qname.namespace else ns

    /*
     * Lookup the tag name in all enums in the current namespace.
     */
    val namespaceMatches = mutable.Set.empty[NamedAst.Declaration.Enum]
    for ((enumName, decl) <- program.enums.getOrElse(namespace, Map.empty[String, NamedAst.Declaration.Enum])) {
      for ((tagName, caze) <- decl.cases) {
        if (tag.name == tagName) {
          namespaceMatches += decl
        }
      }
    }

    // Case 2.1: Exact match found in namespace. Simply return it.
    if (namespaceMatches.size == 1) {
      return Ok(namespaceMatches.head)
    }

    // Case 2.2: No matches found in namespace.
    if (namespaceMatches.isEmpty) {
      return Err(TypeError.UnresolvedTag(qname, tag, ns, tag.loc))
    }

    // Case 2.3: Multiple matches found in namespace...
    ???
  }

}
