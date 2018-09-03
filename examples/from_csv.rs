extern crate csv;
extern crate dcmimu;

use dcmimu::DCMIMU;
use std::io;
use std::str::FromStr;

fn parse_float(s: &str) -> f64 {
    f64::from_str(s).unwrap()
}

fn parse_int(s: &str) -> u32 {
    u32::from_str(s).unwrap()
}

fn time(sr: &csv::StringRecord) -> f64 {
    let tv_sec = parse_int(sr.get(0).unwrap());
    let tv_usec = parse_int(sr.get(1).unwrap());
    tv_sec as f64 + tv_usec as f64 / 1000000.0
}

fn acc(sr: &csv::StringRecord) -> (f64, f64, f64) {
    (parse_float(sr.get(7).unwrap()),
    parse_float(sr.get(8).unwrap()),
    parse_float(sr.get(9).unwrap()))
}

fn gyro(sr: &csv::StringRecord) -> (f64, f64, f64) {
    (parse_float(sr.get(4).unwrap()),
    parse_float(sr.get(5).unwrap()),
    parse_float(sr.get(6).unwrap()))
}

fn reference(sr: &csv::StringRecord) -> (f64, f64, f64) {
    (parse_float(sr.get(11).unwrap()),
    parse_float(sr.get(12).unwrap()),
    parse_float(sr.get(13).unwrap()))
}

// CSV should be formated as follows:
// "tv_sec","tv_usec","data_time","data_index","w_x","w_y","w_z",
// "acc_x","acc_y","acc_z","temperature","pitch","yaw","roll",
fn main() {
    let mut dcmimu = DCMIMU::new();
    let mut rdr = csv::Reader::from_reader(io::stdin());
    let mut prev_t: f64 = 0.0;
    for result in rdr.records() {
        let record = result.unwrap();
        let time = time(&record);
        let (ax, ay, az) = acc(&record);
        let (gx, gy, gz) = gyro(&record);
        let (p, y, r) = reference(&record);
        let dt = if prev_t == 0.0 {
            0.00
        } else {
            time - prev_t
        };
        prev_t = time;
        dcmimu.update((gx as f32, gy as f32, gz as f32),
            (ax as f32, ay as f32, az as f32), dt as f32);
        println!("{:2.8} {:2.8} {:2.8} | {:2.8} {:2.8} {:2.8}",
            dcmimu.yaw(), dcmimu.pitch(), dcmimu.roll(),
            y, p, r);
    }
}

