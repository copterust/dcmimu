extern crate csv;
extern crate dcmimu;

use dcmimu::DCMIMU;
use std::io;
use std::str::FromStr;

fn parse_float(s: &str) -> f32 {
    f32::from_str(s).unwrap()
}

fn parse_int(s: &str) -> u32 {
    u32::from_str(s).unwrap()
}

fn time(sr: &csv::StringRecord) -> f64 {
    let tv_sec = parse_int(sr.get(0).unwrap());
    let tv_usec = parse_int(sr.get(1).unwrap());
    tv_sec as f64 + tv_usec as f64 / 1000000.0
}

fn acc(sr: &csv::StringRecord) -> (f32, f32, f32) {
    (
        parse_float(sr.get(5).unwrap()),
        parse_float(sr.get(6).unwrap()),
        parse_float(sr.get(7).unwrap()),
    )
}

fn gyro(sr: &csv::StringRecord) -> (f32, f32, f32) {
    (
        parse_float(sr.get(2).unwrap()),
        parse_float(sr.get(3).unwrap()),
        parse_float(sr.get(4).unwrap()),
    )
}

fn reference(sr: &csv::StringRecord) -> (f32, f32, f32) {
    (
        parse_float(sr.get(8).unwrap()),
        parse_float(sr.get(9).unwrap()),
        parse_float(sr.get(10).unwrap()),
    )
}

// CSV should be formated as follows:
// "tsecs","tnano","gx","gy","gz","ax","ay","az","yaw","pitch","roll",
// where yaw,pitch,roll are expected results to compare againts.
// Example:
// https://gist.github.com/b979c815a5b4feddbff16185cdf57bcf
fn main() {
    let mut dcmimu = DCMIMU::new();
    let mut rdr = csv::Reader::from_reader(io::stdin());
    let mut prev_t: f64 = 0.0;
    for result in rdr.records() {
        let record = result.unwrap();
        let time = time(&record);
        let (ax, ay, az) = acc(&record);
        let (gx, gy, gz) = gyro(&record);
        let (y, p, r) = reference(&record);
        let dt = if prev_t == 0.0 { 0.00 } else { time - prev_t };
        prev_t = time;
        dcmimu.update((gx, gy, gz), (ax, ay, az), dt as f32);
        let euer_angles = dcmimu.to_euler_angles();
        println!(
            "{:2.8},{:2.8},{:2.8},{:2.8},{:2.8},{:2.8}",
            euer_angles.yaw,
            euer_angles.pitch,
            euer_angles.roll,
            y,
            p,
            r
        );
    }
}
