/*
 * Copyright 2017 Liam Palmer, Esben Bjerre
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

use RedBlackTree.RedBlackTree;

///
/// The Set type.
///
/// A set is currently represented internally as a red-black tree.
///
enum Set[t] {
    case Set(RedBlackTree[t, Unit])
}

instance Boxable[Set[t]] with Order[t], ToString[t] {
    pub def box(x: Set[t]): Boxed = {
        let value = x as ##java.lang.Object;
        let compare = (o1, o2) -> Order.compare(o1 as Set[t], o2 as Set[t]);
        let toString = o -> ToString.toString(o as Set[t]);
        BoxedObject(value, compare, toString)
    }
    pub def unbox(x: Boxed): Set[t] = match x {
        case BoxedObject(v, _, _) => v as Set[t]
        case _                    => ?bug
    }
}

instance ToString[Set[a]] with ToString[a] {
    pub def toString(xs: Set[a]): String =
        let s = Set.foldRight((x, acc) -> "${x}, ${acc}", "", xs);
        "Set#{${String.dropRight(2, s)}}"
}

instance Hash[Set[a]] with Hash[a] {
    pub def hash(xs: Set[a]): Int =
        Set.foldLeft((acc, x) -> acc * 37 + Hash.hash(x), 7793, xs)
}

instance Eq[Set[a]] with Eq[a], Order[a] {
    pub def eq(xs: Set[a], ys: Set[a]): Bool =
        Set.isSubsetOf(xs, ys) and Set.isSubsetOf(ys, xs)
}

instance Order[Set[a]] with Order[a] {

    pub def compare(xs: Set[a], ys: Set[a]): Comparison =
        Set.toList(xs) <=> Set.toList(ys)

}

instance Foldable[Set] {
    pub def foldLeft(f: (b, a) -> b & e, s: b, xs: Set[a]): b & e = Set.foldLeft(f, s, xs)
    pub def foldRight(f: (a, b) -> b & e, s: b, xs: Set[a]): b & e = Set.foldRight(f, s, xs)
}

namespace Set {

    ///
    /// Returns the size of `xs`.
    ///
    @Time(size(xs)) @Space(1)
    pub def size(xs: Set[a]): Int32 =
        let Set(es) = xs;
        RedBlackTree.size(es)






    ///
    /// Returns the empty set.
    ///
    /// `Set#{}` is syntactic sugar for `empty` (`Set#{} = empty()`).
    ///
    @Time(1) @Space(1)
    pub def empty(): Set[a] = Set(RedBlackTree.empty())

    ///
    /// Returns the singleton set containing `x`.
    ///
    /// `Set#{x}` is syntactic sugar for `singleton` (`Set#{x} = singleton(x)`).
    ///
    @Time(1) @Space(1)
    pub def singleton(x: a): Set[a] with Order[a] = insert(x, empty())

    ///
    /// Returns a set of all integers between `b` (inclusive) and `e` (exclusive).
    ///
    /// Returns `empty()` if `b >= e`.
    ///
    @Time(e - b) @Space(e - b)
    pub def range(b: Int32, e: Int32): Set[Int32] =
        if (b >= e) empty() else rangeHelper(b, e - 1, empty())

    ///
    /// Adds `x` to `xs`.
    ///
    @Time(Int32.log2(size(xs))) @Space(Int32.log2(size(xs)))
    pub def insert(x: a, xs: Set[a]): Set[a] with Order[a] =
        let Set(es) = xs;
        Set(RedBlackTree.insert(x, (), es))

    ///
    /// Removes `x` from `xs`.
    ///
    @Time(Int32.log2(size(xs))) @Space(Int32.log2(size(xs)))
    pub def delete(x: a, xs: Set[a]): Set[a] with Order[a] =
        let Set(es) = xs;
        Set(RedBlackTree.delete(x, es))





    ///
    /// Returns true if and only if `xs` is the empty set.
    ///
    @Time(1) @Space(1)
    pub def isEmpty(xs: Set[a]): Bool =
        let Set(es) = xs;
        RedBlackTree.isEmpty(es)

    ///
    /// Returns true if and only if `x` is a member of `xs`.
    ///
    @Time(Int32.log2(size(xs))) @Space(1)
    pub def memberOf(x: a, xs: Set[a]): Bool with Order[a] =
        let Set(es) = xs;
        RedBlackTree.memberOf(x, es)

    ///
    /// Optionally finds the smallest element of `xs` according to the `Order` on `a`.
    ///
    /// Returns `None` if `xs` is empty.
    ///
    pub def minimum(xs: Set[a]): Option[a] =
        let Set(es) = xs;
        RedBlackTree.minimum(es) |> Option.map(fst)

    ///
    /// Optionally finds the smallest element of `xs` according to the given comparator `cmp`.
    ///
    /// Returns `None` if `xs` is empty.
    ///
    pub def minimumBy(cmp: (a, a) -> Comparison, xs: Set[a]): Option[a] =
        reduceLeft(Order.minBy(cmp), xs)

    ///
    /// Optionally finds the largest element of `xs` according to the `Order` on `a`.
    ///
    /// Returns `None` if `xs` is empty.
    ///
    pub def maximum(xs: Set[a]): Option[a] =
        let Set(es) = xs;
        RedBlackTree.maximum(es) |> Option.map(fst)

    ///
    /// Optionally finds the largest element of `xs` according to the given comparator `cmp`.
    ///
    /// Returns `None` if `xs` is empty.
    ///
    pub def maximumBy(cmp: (a, a) -> Comparison, xs: Set[a]): Option[a] =
        reduceLeft(Order.maxBy(cmp), xs)

    ///
    /// Returns true if and only if every element in `xs` appears in `ys`.
    ///
    @Time(size(xs)) @Space(size(xs) * Int32.log2(size(ys)))
    pub def isSubsetOf(xs: Set[a], ys: Set[a]): Bool with Order[a] = forall(x -> memberOf(x, ys), xs)

    ///
    /// Returns true if and only if every element in `xs` appears in `ys` and `xs != ys`.
    ///
    @Time(size(xs)) @Space(size(xs) * Int32.log2(size(ys)))
    pub def isProperSubsetOf(xs: Set[a], ys: Set[a]): Bool with Order[a] = size(xs) != size(ys) and isSubsetOf(xs, ys)

    ///
    /// Alias for `findLeft`.
    ///
    /// The function `f` must be pure.
    ///
    @Time(time(f) * size(xs)) @Space(space(f) * Int32.log2(size(xs)))
    pub def find(f: a -> Bool, xs: Set[a]): Option[a] = findLeft(f, xs)

    ///
    /// Optionally returns the first element of `xs` that satisfies the predicate `f` when searching from left to right.
    ///
    /// The function `f` must be pure.
    ///
    @Time(time(f) * size(xs)) @Space(space(f) * Int32.log2(size(xs)))
    pub def findLeft(f: a -> Bool, xs: Set[a]): Option[a] =
        let Set(es) = xs;
        RedBlackTree.findLeft((x, _) -> f(x), es) |> Option.map(fst)

    ///
    /// Optionally returns the first element of `xs` that satisfies the predicate `f` when searching from right to left.
    ///
    /// The function `f` must be pure.
    ///
    @Time(time(f) * size(xs)) @Space(space(f) * Int32.log2(size(xs)))
    pub def findRight(f: a -> Bool, xs: Set[a]): Option[a] =
        let Set(es) = xs;
        RedBlackTree.findRight((x, _) -> f(x), es) |> Option.map(fst)






    ///
    /// Returns the result of applying `combine` to all the elements in `xs`, using `empty` as the initial value.
    ///
    pub def fold(xs: Set[a]): a with Monoid[a] = Foldable.fold(xs)

    ///
    /// Applies `f` to a start value `s` and all elements in `xs` going from left to right.
    ///
    /// That is, the result is of the form: `f(...f(f(s, x1), x2)..., xn)`.
    ///
    pub def foldLeft(f: (b, a) -> b & e, s: b, xs: Set[a]): b & e =
        let Set(es) = xs;
        RedBlackTree.foldLeft((acc, k, _) -> f(acc, k), s, es)

    ///
    /// Applies `f` to a start value `s` and all elements in `xs` going from right to left.
    ///
    /// That is, the result is of the form: `f(x1, ...f(xn-1, f(xn, s))...)`.
    ///
    pub def foldRight(f: (a, b) -> b & e, s: b, xs: Set[a]): b & e =
        let Set(es) = xs;
        RedBlackTree.foldRight((k, _, acc) -> f(k, acc), s, es)

    ///
    /// Applies `f` to all elements in `xs` going from left to right until a single value `v` is obtained.  Returns `Some(v)`.
    /// That is, the result is of the form: `Some(f(...f(f(x1, x2), x3)..., xn))`
    /// Returns `None` if `xs` is the empty set.
    ///
    pub def reduceLeft(f: (a, a) -> a & e, xs: Set[a]): Option[a] & e =
        let Set(es) = xs;
        RedBlackTree.reduceLeft((x, _, y, _) -> (f(x, y), ()), es) |> Option.map(fst)

    ///
    /// Applies `f` to all elements in `xs` going from right to left until a single value `v` is obtained.  Returns `Some(v)`.
    /// That is, the result is of the form: `Some(f(x1, ...f(xn-2, f(xn-1, xn))...))`
    /// Returns `None` if `xs` is the empty set.
    ///
    pub def reduceRight(f: (a, a) -> a & e, xs: Set[a]): Option[a] & e =
        let Set(es) = xs;
        RedBlackTree.reduceRight((x, _, y, _) -> (f(x, y), ()), es) |> Option.map(fst)






    ///
    /// Returns the number of elements in `xs` that satisfy the predicate `f`.
    ///
    /// The function `f` must be pure.
    ///
    @Time(time(f) * size(xs)) @Space(space(f) * Int32.log2(size(xs)))
    pub def count(f: a -> Bool, xs: Set[a]): Int32 =
        foldLeft((acc, x) -> if (f(x)) 1 + acc else acc, 0, xs)

    ///
    /// Returns the union of the elements in `xs`.
    ///
    @Time(size(xs) * Int32.log2(size(xs))) @Space(size(xs) * Int32.log2(size(xs)))
    pub def flatten(xs: Set[Set[a]]): Set[a] with Order[a] =
        foldLeft((acc, x) -> union(acc, x), empty(), xs)

    ///
    /// Returns `true` if and only if at least one element in `xs` satisfies the predicate `f`.
    ///
    /// Returns `false` if `xs` is the empty set.
    ///
    /// The function `f` must be pure.
    ///
    @Time(time(f) * size(xs)) @Space(space(f) * Int32.log2(size(xs)))
    pub def exists(f: a -> Bool, xs: Set[a]): Bool =
        let Set(es) = xs;
        RedBlackTree.exists((x, _) -> f(x), es)

    ///
    /// Returns `true` if and only if all elements in `xs` satisfy the predicate `f`.
    ///
    /// Returns `true` if `xs` is the empty set.
    ///
    /// The function `f` must be pure.
    ///
    @Time(time(f) * size(xs)) @Space(space(f) * Int32.log2(size(xs)))
    pub def forall(f: a -> Bool, xs: Set[a]): Bool =
        let Set(es) = xs;
        RedBlackTree.forall((x, _) -> f(x), es)






    ///
    /// Returns the union of `xs` and `ys`.
    ///
    @Time(Int32.min(size(xs), size(ys)) * Int32.log2(size(xs) + size(ys)))
    @Space(Int32.min(size(xs), size(ys)) * Int32.log2(size(xs) + size(ys)))
    pub def union(xs: Set[a], ys: Set[a]): Set[a] with Order[a] =
        use RedBlackTree.{blackHeight, foldLeft, insert};
        let Set(as) = xs;
        let Set(bs) = ys;
        if (blackHeight(as) < blackHeight(bs))
            Set(foldLeft((acc, x, _) -> insert(x, (), acc), bs, as))
        else
            Set(foldLeft((acc, y, _) -> insert(y, (), acc), as, bs))

    ///
    /// Returns the intersection of `xs` and `ys`.
    ///
    @Time(size(ys) * Int32.log2(size(xs))) @Space(size(ys) * Int32.log2(size(xs)))
    pub def intersection(xs: Set[a], ys: Set[a]): Set[a] with Order[a] =
        foldLeft((acc, x) -> if (memberOf(x, xs)) insert(x, acc) else acc, empty(), ys)

    ///
    /// Returns the difference of `xs` and `ys`, i.e. `xs - ys`.
    ///
    @Time(size(xs) * Int32.log2(size(ys))) @Space(size(xs) * Int32.log2(size(ys)))
    pub def difference(xs: Set[a], ys: Set[a]): Set[a] with Order[a] =
        foldLeft((acc, x) -> if (not memberOf(x, ys)) insert(x, acc) else acc, empty(), xs)

    ///
    /// Returns all subsets of `xs`.
    ///
    @Time(size(xs) * (1 <<< size(xs))) @Space(size(xs) * (1 <<< size(xs)))
    pub def subsets(xs: Set[a]): Set[Set[a]] with Order[a] =
        foldLeft((acc, x) -> union(map(y -> insert(x, y), acc), acc), insert(empty(), empty()), xs)






    ///
    /// Returns the set of all elements of `xs` that satisfy the predicate `f`.
    ///
    /// The function `f` must be pure.
    ///
    @Time(time(f) * size(xs)) @Space(space(f) * size(xs))
    pub def filter(f: a -> Bool, xs: Set[a]): Set[a] with Order[a] =
        foldLeft((acc, x) -> if (f(x)) insert(x, acc) else acc, empty(), xs)

    ///
    /// Returns the result of applying `f` to every element in `xs`.
    ///
    /// Note: The returned set may be smaller than `xs`.
    ///
    @Time(time(f) * size(xs)) @Space(space(f) * size(xs))
    pub def map(f: a -> b & e, xs: Set[a]): Set[b] & e with Order[b] =
        foldLeft((acc, x) -> insert(f(x), acc), empty(), xs)

    ///
    /// Returns the result of applying `f` to every element in `xs` and taking the union.
    ///
    @Time(time(f) * size(xs) * Int32.log2(size(xs))) @Space(space(f) * size(xs) * Int32.log2(size(xs)))
    pub def flatMap(f: a -> Set[b] & e, xs: Set[a]): Set[b] & e with Order[b] =
        foldLeft((acc, x) -> union(acc, f(x)), empty(), xs)

    ///
    /// Replaces the element `x` with `y` if `x` is in `xs`. Otherwise, returns `xs`.
    ///
    /// Note: The returned set may be smaller than `xs`.
    ///
    @Time(Int32.log2(size(xs))) @Space(Int32.log2(size(xs)))
    pub def replace(x: a, y: a, xs: Set[a]): Set[a] with Order[a] =
        if (memberOf(x, xs)) insert(y, delete(x, xs)) else xs

    ///
    /// Returns a pair of sets `(ys, zs)`.
    ///
    /// `ys` contains all elements of `xs` that satisfy the predicate `f`.
    /// `zs` contains all elements of `xs` that do not satisfy the predicate `f`.
    ///
    /// The function `f` must be pure.
    ///
    @Time(size(xs) * Int32.log2(size(xs))) @Space(size(xs) * Int32.log2(size(xs)))
    pub def partition(f: a -> Bool, xs: Set[a]): (Set[a], Set[a]) with Order[a] =
        foldLeft((acc, x) ->
            let (a, b) = acc;
            if (f(x))
                (insert(x, a), b)
            else
                (a, insert(x, b)), (empty(), empty()),
            xs)




    ///
    /// Returns `xs` as a mutable set.
    ///
    @Time(1) @Space(1)
    pub def toMut(xs: Set[a]): MutSet[a] & Impure =
        MutSet(ref xs)

    ///
    /// Returns the set `xs` as a list.
    ///
    @Time(size(xs)) @Space(size(xs))
    pub def toList(xs: Set[a]): List[a] =
        foldLeft((acc, x) -> x :: acc, Nil, xs)

    ///
    /// Returns the association set `xs` as a map.
    ///
    /// If `xs` contains multiple mappings with the same key, `toMap` does not
    /// make any guarantees about which mapping will be in the resulting map.
    ///
    @Time(size(xs)) @Space(size(xs))
    pub def toMap(xs: Set[(a, b)]): Map[a, b] with Order[a] =
        foldRight((x, acc) -> let (k, v) = x; Map.insert(k, v, acc), Map#{}, xs)

    ///
    /// Applies `f` to every element of `xs`.
    ///
    /// The function `f` must be impure.
    ///
    @Time(size(f) * size(xs)) @Space(space(f) * Int32.log2(size(xs)))
    pub def foreach(f: a ~> Unit, xs: Set[a]): Unit & Impure =
        let Set(es) = xs;
        RedBlackTree.foreach((x, _) -> f(x), es)

    ///
    /// Build a set by applying `f` to the seed value `st`.
    ///
    /// `f` should return `Some(a,st1)` to signal a new set element `a` and a new seed value `st1`.
    ///
    /// `f` should return `None` to signal the end of building the set.
    ///
    pub def unfold(f: s -> Option[(a, s)] & e, st: s): Set[a] & e with Order[a] = unfoldHelper(f, st, empty())

    ///
    /// Helper function for `unfold`.
    ///
    def unfoldHelper(f: s -> Option[(a, s)] & e, st: s, ac: Set[a]): Set[a] & e with Order[a] = match f(st) {
        case None => ac
        case Some((a, st1)) => unfoldHelper(f, st1, insert(a, ac))
    }

    ///
    /// Build a set by applying the function `next` to `()`. `next` is expected to encapsulate
    /// a stateful resource such as a file handle that can be iterated.
    ///
    /// `next` should return `Some(a)` to signal a value pair `a`.
    ///
    /// `next` should return `None` to signal the end of building the set.
    ///
    pub def unfoldWithIter(next: () ~> Option[a]): Set[a] & Impure with Order[a] = unfoldWithIterHelper(next, empty())

    ///
    /// Helper function for `unfoldWithIter`.
    ///
    def unfoldWithIterHelper(next: () ~> Option[a], ac: Set[a]): Set[a] & Impure with Order[a] = match next() {
        case None => ac
        case Some(a) => unfoldWithIterHelper(next, insert(a, ac))
    }

    ///
    /// Extracts a range of elements from the set `xs`.
    ///
    /// That is, the result is a list of all elements `x` where `p(x)` returns `EqualTo`.
    ///
    pub def query(p: a -> Comparison & e, xs: Set[a]): List[a] & e =
        let Set(es) = xs;
        RedBlackTree.query(p, (k, _) -> k, es)

    ///
    /// Applies `f` to all elements from the set `xs` where `p(x)` returns `EqualTo`.
    ///
    /// The function `f` must be impure.
    ///
    pub def queryWith(p: a -> Comparison & e, f: a ~> Unit, xs: Set[a]): Unit & Impure =
        let Set(es) = xs;
        RedBlackTree.queryWith(p, (k, _) -> f(k), es)

    ///
    /// Helper function for `range`.
    ///
    def rangeHelper(b: Int32, n: Int32, acc: Set[Int32]): Set[Int32] =
        if (n < b) acc else rangeHelper(b, n - 1, insert(n, acc))

}
