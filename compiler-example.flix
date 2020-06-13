Interpreter("x86").
Compiler("C++", "C", "x86").
Compiler("C", "x86", "x86").
Compiler("Java", "x86", "C++").

Compiler(a, c, lang) :-
    Compiler(a, b, lang),
    Compiler(b, c, lang),
    Interpreter(lang).

Compiler(src, dst, lang2) :-
    Compiler(src, dst, lang1),
    Compiler(lang1, lang2, lang3),
    Interpreter(lang3).