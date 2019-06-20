use cgmath;
use std::f64::consts::PI;
use std::f64::EPSILON;
use std::f64;
use libm::ldexp;
use std::ops::Mul;
use std::ops::Add;

use rand;
use rand::Rng;
use s2::cap::Cap;
use s2::cellid::*;
use s2::point::{self, Point};
use r3;

pub static DBL_EPSILON: f64 = 2.2204460492503131e-16;

pub fn rng() -> rand::prng::XorShiftRng {
    use rand::prelude::*;

    rand::prng::XorShiftRng::from_rng(rand::thread_rng()).expect("failed to get rng")
}

/// skewed_int returns a number in the range [0,2^max_log-1] with bias towards smaller numbers.
pub fn skewed_int<R>(rng: &mut R, max_log: usize) -> usize
where
    R: Rng,
{
    let base = rng.gen_range(0, max_log + 1);
    rng.gen_range(0, 1 << 31) & ((1 << base) - 1)
}

/// cap returns a cap with a random axis such that the log of its area is
/// uniformly distributed between the logs of the two given values. The log of
/// the cap angle is also approximately uniformly distributed.
pub fn cap<R>(rng: &mut R, min_area: f64, max_area: f64) -> Cap
where
    R: Rng,
{
    let cap_area = max_area * (min_area / max_area).powf(rng.gen_range(0., 1.));
    Cap::from_center_area(&point(rng), cap_area)
}

/// point returns a random unit-length vector.
pub fn point<R: Rng>(rng: &mut R) -> Point {
    Point::from_coords(
        rng.gen_range(-1., 1.),
        rng.gen_range(-1., 1.),
        rng.gen_range(-1., 1.),
    )
}

pub fn frame<R: Rng>(rng: &mut R) -> cgmath::Matrix3<f64> {
    let z = point(rng);
    frame_at_point(rng, z)
}

pub fn frame_at_point<R: Rng>(rng: &mut R, z: Point) -> cgmath::Matrix3<f64> {
    let p = point(rng);
    let x = z.cross(&p).normalize();
    let y = z.cross(&x).normalize();

    cgmath::Matrix3::from_cols(x.into(), y.into(), z.into())
}

pub fn cellid<R>(rng: &mut R) -> CellID
where
    R: Rng,
{
    let level = rng.gen_range(0, MAX_LEVEL + 1);
    cellid_for_level(rng, level)
}

pub fn cellid_for_level<R>(rng: &mut R, level: u64) -> CellID
where
    R: Rng,
{
    let face = rng.gen_range(0, NUM_FACES as u64);
    let pos = rng.next_u64() & ((1 << POS_BITS) - 1);
    let cellid = CellID::from_face_pos_level(face, pos, level);
    assert_eq!(face, cellid.face() as u64);
    assert_eq!(level, cellid.level());

    cellid
}

pub fn one_in<R>(rng: &mut R, n: u64) -> bool
where
    R: Rng,
{
    rng.gen_range(0, n) == 0
}

// sample_point_from_cap returns a point chosen uniformly at random (with respect
// to area) from the given cap.
pub fn sample_point_from_cap<R>(rng: &mut R, c: Cap) -> Point
where
    R: Rng,
{
    // We consider the cap axis to be the "z" axis. We choose two other axes to
    // complete the coordinate frame.
    let center = c.center();
    let m = center.frame();

    // The surface area of a spherical cap is directly proportional to its
    // height. First we choose a random height, and then we choose a random
    // point along the circle at that height.
    let h = rng.gen_range(0., 1.) * c.height();
    let theta = 2. * PI * rng.gen_range(0., 1.);
    let r = (h * (2. - h)).sqrt();

    // The result should already be very close to unit-length, but we might as
    // well make it accurate as possible.
    point::from_frame(
        &m,
        &Point::from_coords(theta.cos() * r, theta.sin() * r, 1. - h).normalize(),
    )
}

// randomBits returns a 64-bit random unsigned integer whose lowest "num" are random, and
// whose other bits are zero.
pub fn random_bits<R>(rng: &mut R, num: u32) -> u64
where
    R: Rng,
{
    // Make sure the request is for not more than 63 bits.
    let mut n = num;
    if n > 63 {
        n = 63
    }
    return rng.gen::<u64>() & ((1 << n) - 1)
}

// OriginPoint returns a unique "origin" on the sphere for operations that need a fixed
// reference point. In particular, this is the "point at infinity" used for
// point-in-polygon testing (by counting the number of edge crossings).
//
// It should *not* be a point that is commonly used in edge tests in order
// to avoid triggering code to handle degenerate cases (this rules out the
// north and south poles). It should also not be on the boundary of any
// low-level S2Cell for the same reason.
pub fn origin_point() -> Point {
    return Point(r3::vector::Vector{x: -0.0099994664350250197, y: 0.0025924542609324121, z: 0.99994664350250195})
}


// PointFromCoords creates a new normalized point from coordinates.
//
// This always returns a valid point. If the given coordinates can not be normalized
// the origin point will be returned.
//
// This behavior is different from the C++ construction of a S2Point from coordinates
// (i.e. S2Point(x, y, z)) in that in C++ they do not Normalize.
pub fn point_from_coords(x:f64, y:f64, z:f64) -> Point {
    if x == 0.0 && y == 0.0 && z == 0.0 {
        return origin_point()
    }
    return Point(r3::vector::Vector{x, y, z}.normalize())
}

// randomPoint returns a random unit-length vector.
pub fn random_point() -> Point {
    return point_from_coords(random_uniform_f64(-1.0, 1.0),
        random_uniform_f64(-1.0, 1.0), random_uniform_f64(-1.0, 1.0))
}

// randomUniformFloat64 returns a uniformly distributed value in the range [min, max).
pub fn random_uniform_f64(min:f64, max:f64) -> f64 {
    return min + random_f64() * (max - min)
}

// randomUniformInt returns a uniformly distributed integer in the range [0,n).
// NOTE: This is replicated here to stay in sync with how the C++ code generates
// uniform randoms. (instead of using Go's math/rand package directly).
pub fn random_uniform_int(n: i64) -> i64 {
    return (random_f64() * n as f64) as i64
}

// randomFloat64 returns a uniformly distributed value in the range [0,1).
// Note that the values returned are all multiples of 2**-53, which means that
// not all possible values in this range are returned.
pub fn random_f64() -> f64 {
    let mut rng = rng();

    const RANDOM_FLOAT_BITS: i32 = 53;
    return ldexp(random_bits(&mut rng, RANDOM_FLOAT_BITS as u32) as f64, -RANDOM_FLOAT_BITS)
}

// float64Near reports whether the two values are within the given epsilon.
pub fn f64_near(x: f64, y: f64, e: f64) -> bool {
    return (x-y).abs() <= e
}

// float64Eq reports whether the two values are within the default epsilon.
pub fn f64_eq(x: f64, y: f64) -> bool { return f64_near(x, y, EPSILON) }

// perturbedCornerOrMidpoint returns a Point from a line segment whose endpoints are
// difficult to handle correctly. Given two adjacent cube vertices P and Q,
// it returns either an edge midpoint, face midpoint, or corner vertex that is
// in the plane of PQ and that has been perturbed slightly. It also sometimes
// returns a random point from anywhere on the sphere.
pub fn perturbed_corner_or_midpoint(p: Point, q: Point) -> Point {
    let mut a = p.mul((random_uniform_int(3) - 1) as f64).add(q.mul((random_uniform_int(3) - 1) as f64));
    let mut rng = rng();
    if one_in(&mut rng, 10) {
        // This perturbation often has no effect except on coordinates that are
        // zero, in which case the perturbed value is so small that operations on
        // it often result in underflow.
        a = a.add(random_point().mul(1e-300_f64.powf(random_f64())));
    } else if one_in(&mut rng, 2) {
        // For coordinates near 1 (say > 0.5), this perturbation yields values
        // that are only a few representable values away from the initial value.
        a = a.add(random_point().mul(4.0 * DBL_EPSILON));
    } else {
        // A perturbation whose magnitude is in the range [1e-25, 1e-10].
        a = a.add(random_point().mul(1e-10 * 1e-15_f64.powf(random_f64())));
    }

    if a.0.norm2() < f64::MIN {
        // If a.Norm2() is denormalized, Normalize() loses too much precision.
        return perturbed_corner_or_midpoint(p, q)
    }
    return a
}
