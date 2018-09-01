# `dcmimu`

> An algorithm for fusing low-cost triaxial MEMS gyroscope and accelerometer measurements.

A `no_std` Rust port of [the original](https://github.com/hhyyti/dcm-imu).

[![Build Status](https://travis-ci.org/copterust/dcmimu.svg?branch=master)](https://travis-ci.org/copterust/dcmimu)

*NOTE*: `libm` still [doesn't with overflow checks](https://github.com/japaric/libm/issues/4),
so you have to compile your project with `--release`.
Leave a comment in the linked issue to [raise awareness](https://www.youtube.com/watch?v=KbZIFZm204E).

## Credentials

[Heikki Hyyti and Arto Visala, "A DCM Based Attitude Estimation Algorithm for Low-Cost MEMS IMUs," International Journal of Navigation and Observation, vol. 2015, Article ID 503814, 18 pages, 2015](http://dx.doi.org/10.1155/2015/503814).

## Usage

Library is available via crates.io [![crates.io](http://meritbadge.herokuapp.com/dcmimu?style=flat-square)](https://crates.io/crates/dcmimu).

```rust

# Create DCMIMU:
let mut dcmimu = DCMIMU::new();
let mut prev_t_ms = now();
loop {
    # get gyroscope and accelerometer measurement from your sensors:
    let gyro = sensor.read_gyro();
    let accel = sensor.read_accel();
    # Convert measurements to SI if needed.
    # Get time difference since last update:
    let t_ms = now();
    let dt_ms = t_ms - prev_t_ms
    prev_t_ms = t_ms
    # Update dcmimu states (don't forget to use SI):
    let dcm = dcmimu.update((gyro.x, gyro.y, gyro.z),
                            (accel.x, accel.y, accel.z),
                            dt_ms.seconds());
    println!("Roll: {}; yaw: {}; pitch: {}", dcm.roll, dcm.yaw, dcm.pitch);
    # Measurements can also be queried without updating:
    println!("{:?} == {}, {}, {}", dcmimu.all(), dcmimu.roll(), dcmimu.yaw(), dcmimu.pitch());
}

```

Check out [mpu9250](https://crates.io/crates/mpu9250) for accelerometer/gyroscrope sensors driver.

## Documentation

Available via [docs.rs](https://docs.rs/dcmimu/).

## License

[MIT license](http://opensource.org/licenses/MIT).
