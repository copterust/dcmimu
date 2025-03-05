extern crate dcmimu;

use dcmimu::DCMIMU;
use std::io;
use std::io::BufRead;
use std::str::FromStr;

fn parse_float(s: &str) -> f32 {
    f32::from_str(s).unwrap()
}

// format:
// ax,ay,az,gx,gy,gz,dt_s,y,p,r
fn main() {
    let mut dcmimu = DCMIMU::new();
    let stdin = io::stdin();
    for line in stdin.lock().lines() {
        let l = match line {
            Ok(ll) => ll,
            Err(_e) => return,
        };
        let split = l.split(" ");
        let vec = split.collect::<Vec<&str>>();
        let (ax, ay, az) = (
            parse_float(vec[0]),
            parse_float(vec[1]),
            parse_float(vec[2]),
        );
        let (gx, gy, gz) = (
            parse_float(vec[3]),
            parse_float(vec[4]),
            parse_float(vec[5]),
        );
        let dt_s = parse_float(vec[6]);
        let (ry, rp, rr) = (
            parse_float(vec[7]),
            parse_float(vec[8]),
            parse_float(vec[9]),
        );
        dcmimu.update((gx, gy, gz), (ax, ay, az), dt_s);
        let ypr = dcmimu.to_euler_angles();

        for f in [ypr.yaw, ypr.pitch, ypr.roll, ry, rp, rr].into_iter() {
            let mut b = ryu::Buffer::new();
            let s = b.format(*f);
            print!("{}, ", s);
        }
        println!();
    }
}
