extern crate dcmimu;

use dcmimu::DCMIMU;

fn main() {
    let mut imu = DCMIMU::new();
    imu.update_only((0.0, 0.0, 0.0), (0.0, 0.0, 0.0), 0.1);
    let dcm = imu.to_euler_angles();
}
