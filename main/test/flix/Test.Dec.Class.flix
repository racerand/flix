namespace Test/Dec/Class {
    namespace Test01 {

        class X[a]
    }

    namespace Test02 {

        lawless class X[a] {
            pub def f(x: a): Bool
        }

    }

    namespace Test03 {
        lawless class X[a] {
            pub def f(x: a): Bool
        }

        instance X[Int] {
            pub def f(x: Int): Bool = x == 0
        }

        pub def g(x: Int): Bool = Test/Dec/Class/Test03/X.f(x)
    }

    namespace Test04 {
        lawless class C[a] {
            pub def f(x: a): Bool
        }

        lawless class D[a] {
            pub def g(x: a): Bool
        }

        pub def h(x: a): Bool with C[a], D[a] = Test/Dec/Class/Test04/C.f(x) and Test/Dec/Class/Test04/D.g(x)
    }

    namespace Test05 {
        use Test/Dec/Class.List;

        class C[a]

        instance C[List[a]] with C[a]
    }

    namespace Test06 {
        use Test/Dec/Class.List;
        use Test/Dec/Class.List.Nil;
        use Test/Dec/Class.Option;

        pub lawless class C[a] {
            pub def f(x: a): Bool
        }

        instance C[Int] {
            pub def f(x: Int): Bool = x == 0
        }

        instance C[List[a]] with C[a] {
            pub def f(x: List[a]): Bool = {
                match x {
                    case Nil => true
                    case _ => false
                }
            }
        }

        pub def g(x: List[Int]): Bool = Test/Dec/Class/Test06/C.f(x)
    }

    namespace Test07 {
        use Test/Dec/Class/Test06.C;

        pub def g(x: a): a with C[a] = x
    }

    namespace Test08 {
        use Test/Dec/Class.Option;
        use Test/Dec/Class.Option.None;
        use Test/Dec/Class.Option.Some;

        lawless class F[m : Type -> Type] {
            pub def map(f: a -> b, x: m[a]): m[b]
        }

        instance F[Option] {
            pub def map(f: a -> b, x: Option[a]): Option[b] = match x {
                case None => None
                case Some(y) => Some(f(y))
            }
        }
    }

    namespace Test09 {
        lawless class Eff[e : Bool] {
            pub def isPure(f: a -> b & e): Bool
        }

        instance Eff[Pure] {
            pub def isPure(_f: a -> b): Bool = true
        }

        instance Eff[Impure] {
            pub def isPure(_f: a ~> b): Bool = false
        }
    }

    namespace Test10 {
        lawless class C[a] {
            pub def f(x: a): Bool
        }

        lawless class D[a] with C[a]

        lawless class E[a] with D[a]

        instance C[Int] {
            pub def f(_x: Int): Bool = true
        }

        instance D[Int]

        instance E[Int]

        pub def g(x: a): Bool with D[a] = Test/Dec/Class/Test10/C.f(x)

        pub def h(x: a): Bool with E[a] = Test/Dec/Class/Test10/C.f(x)
    }

    namespace Test11 {
        class C[a] {
            pub def f(x: a, y: a): Bool

            // TODO handle namespaces better
            law reflexivity: forall (x: a, y: a) . Test/Dec/Class/Test11/C.f(x, y) == Test/Dec/Class/Test11/C.f(y, x)
        }

        instance C[Int] {
            pub def f(_x: Int, _y: Int): Bool = true
        }
    }

    namespace Test12 {
        lawless class C[a] {
            pub def f(): a

            law l: forall (_x: a) . true

            pub def g(): a

            law m: forall[a: Type](_x: a) . true
        }
    }

    namespace Test13 {
        lawless class C[a] {
            pub def f(): a

            pub def g(): a = Test/Dec/Class/Test13/C.f()

            pub def h(): a = Test/Dec/Class/Test13/C.g()
        }

        instance C[Int] {
            pub def f(): Int = 0

            override pub def h(): Int = 1
        }

        @test
        def testSigOverride01(): Bool = Test/Dec/Class/Test13/C.f() == 0

        @test
        def testSigOverride02(): Bool = Test/Dec/Class/Test13/C.g() == 0

        @test
        def testSigOverride03(): Bool = Test/Dec/Class/Test13/C.h() == 1
    }

    namespace Test14 {
        lawless class C[a] {
            pub def c(x: a): Int
        }

        lawless class D[a] with C[a] {
            pub def d(x: a): Int
        }

        pub def f(x: a): Int with D[a] = Test/Dec/Class/Test14/C.c(x) + Test/Dec/Class/Test14/D.d(x)
    }

    namespace Test15 {
        lawless class C[a]

        lawless class D[a]

        pub def f(x: a, y: b): Bool with C[a], D[b] = ???
    }

    pub enum List[a] {
        case Nil
        case Cons(a, List[a])
    }

    pub enum Option[a] {
        case None
        case Some(a)
    }

    @test
    def empty(): Bool = true
}
