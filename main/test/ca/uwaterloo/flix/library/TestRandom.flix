@test
def testChoose01(): Option[Int32] & Impure =
    let r = Random.new();
    let a = [42i32; 10];
    Random.choose(r, a)

@test
def testChoose02(): Option[String] & Impure =
    let r = Random.new();
    let a = ["Hello World"; 10];
    Random.choose(r, a)

@test
def testChoose03(): Option[Result[Int, Int]] & Impure =
    let r = Random.new();
    let a = [Ok(123), Ok(456), Err(789)];
    Random.choose(r, a)

@test
def testNextBool(): Bool & Impure =
    let r = Random.newWithSeed(0i64);
    Random.nextBool(r)

@test
def testNextFloat32(): Float32 & Impure =
    let r = Random.new();
    Random.nextFloat32(r)

@test
def testNextFloat64(): Float64 & Impure =
    let r = Random.new();
    Random.nextFloat64(r)

@test
def testNextGaussian(): Float64 & Impure =
    let r = Random.new();
    Random.nextGaussian(r)

@test
def testNextInt32(): Int32 & Impure =
    let r = Random.new();
    Random.nextInt32(r)

@test
def testNextInt64(): Int64 & Impure =
    let r = Random.new();
    Random.nextInt64(r)

@test
def testNextNatWithMax(): Int32 & Impure =
    let r = Random.new();
    Random.nextNatWithMax(r, 100)
