/*
 * Copyright 2021 Benjamin Dahse
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

use Fixpoint/Ast.Denotation;

namespace Fixpoint/Ram {
    pub enum RamSym[v] {
        case Full(String, Denotation[v])
        case Delta(String, Denotation[v])
        case New(String, Denotation[v])
    }

    pub def toDenotation(ramSym: RamSym[v]): Denotation[v] = match ramSym {
        case Full(_, den) => den
        case Delta(_, den) => den
        case New(_, den) => den
    }

    instance ToString[RamSym[v]] {
        pub def toString(ramSym: RamSym[v]): String = match ramSym {
            case Full(name, _) => name
            case Delta(name, _) => "Δ${name}"
            case New(name, _) => "Δ${name}'"
        }
    }

    instance Order[RamSym[v]] {
        pub def compare(a: RamSym[v], b: RamSym[v]): Comparison = match a {
            case Full(s1, _) => match b {
                case Full(s2, _) => s1 <=> s2
                case _ => GreaterThan
            }
            case Delta(s1, _) => match b {
                case Delta(s2, _) => s1 <=> s2
                case New(_) => GreaterThan
                case Full(_) => LessThan
            }
            case New(s1, _) => match b {
                case New(s2, _) => s1 <=> s2
                case _ => LessThan
            }
        }
    }
}