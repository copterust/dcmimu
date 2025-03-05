extern crate dcmimu;

use dcmimu::DCMIMU;

fn main() {
    let mut imu = DCMIMU::new();
    imu.update((0.0, 0.0, 0.0), (0.0, 0.0, 0.0), 0.1);
}
