rel Train(src: City, dst: City)
rel Bus(src: City, dst: City)
rel Path(src: City, dst: City)

type alias City = Int

def travel(tg: #{Train| r}, bg: #{Bus| r}, src:City ,dst:City): Option[#{Path}] =
    let p = #{
        Path(x,y) :- Train(x,y).
        Path(x,z) :- Path(x,y), Train(y,z).
        Path(x,y) :- Bus(x,y),not Train(x,y).
        Path(x,z) :- Path(x,y), Bus(y,z),not Train(x,y).
    };
    let m = solve (tg <+> bg <+> p);
    if (m |= Path(src, dst).)
        Some(project Path m)
    else
        None