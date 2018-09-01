//! no_std implementation of Heikki Hyyti and Arto Visala DCM-IMU algorithm.
//! """The DCM-IMU algorithm is designed for fusing low-cost triaxial MEMS
//! gyroscope and accelerometer measurements. An extended Kalman filter is used
//! to estimate attitude in direction cosine matrix (DCM) formation and
//! gyroscope biases online. A variable measurement covariance method is
//! implemented for acceleration measurements to ensure robustness against
//! transient non-gravitational accelerations which usually induce errors
//! to attitude estimate in ordinary IMU-algorithms."""
//!
//! # Usage
//! ```
//! # Create DCMIMU:
//! let mut dcmimu = DCMIMU::new();
//! let mut prev_t_ms = now();
//! loop {
//!     # get gyroscope and accelerometer measurement from your sensors:
//!     let gyro = sensor.read_gyro();
//!     let accel = sensor.read_accel();
//!     # Convert measurements to SI if needed.
//!     # Get time difference since last update:
//!     let t_ms = now();
//!     let dt_ms = t_ms - prev_t_ms
//!     prev_t_ms = t_ms
//!     # Update dcmimu states (don't forget to use SI):
//!     let dcm = dcmimu.update((gyro.x, gyro.y, gyro.z),
//!                             (accel.x, accel.y, accel.z),
//!                             dt_ms.seconds());
//!     println!("Roll: {}; yaw: {}; pitch: {}", dcm.roll, dcm.yaw, dcm.pitch);
//!     # Measurements can also be queried without updating:
//!     println!("{:?} == {}, {}, {}", dcmimu.all(), dcmimu.roll(), dcmimu.yaw(), dcmimu.pitch());
//! }
//! ```
//!

#![no_std]
#![allow(non_snake_case)]
#![deny(warnings)]

extern crate libm;
use libm::{asinf, atan2f, cosf, sinf, sqrtf};

#[cfg_attr(rustfmt, rustfmt_skip)]
pub struct DCMIMU {
    g0: f32,
    g0_2: f32,
    x0: f32, x1: f32, x2: f32, x3: f32, x4: f32, x5: f32,
    q_dcm2: f32,
    q_gyro_bias2: f32,
    r_acc2: f32,
    r_a2: f32,
    a0: f32, a1: f32, a2: f32,
    yaw: f32, pitch: f32, roll: f32,
    // matrix
    P00: f32, P01: f32, P02: f32, P03: f32, P04: f32, P05: f32,
    P10: f32, P11: f32, P12: f32, P13: f32, P14: f32, P15: f32,
    P20: f32, P21: f32, P22: f32, P23: f32, P24: f32, P25: f32,
    P30: f32, P31: f32, P32: f32, P33: f32, P34: f32, P35: f32,
    P40: f32, P41: f32, P42: f32, P43: f32, P44: f32, P45: f32,
    P50: f32, P51: f32, P52: f32, P53: f32, P54: f32, P55: f32,
}

pub const GRAVITY: f32 = 9.81;
const INITIAL_DCM_VARIANCE: f32 = 1.0;
const INITIAL_BIAS_VARIANCE: f32 = 0.1 * 0.1;

impl DCMIMU {
    /// New
    pub fn new() -> Self {
        #[cfg_attr(rustfmt, rustfmt_skip)]
        DCMIMU {
            g0: GRAVITY,
            g0_2: GRAVITY * GRAVITY,
            x0: 0.0,
            x1: 0.0,
            x2: 1.0,
            x3: 0.0,
            x4: 0.0,
            x5: 0.0,
            q_dcm2: 0.1 * 0.1,
            q_gyro_bias2: 0.0001 * 0.0001,
            r_acc2: 0.5 * 0.5,
            r_a2: 10.0 * 10.0,
            a0: 0.0,
            a1: 0.0,
            a2: 0.0,
            yaw: 0.0,
            pitch: 0.0,
            roll: 0.0,
            P00: INITIAL_DCM_VARIANCE,
            P01: 0.0,
            P02: 0.0,
            P03: 0.0,
            P04: 0.0,
            P05: 0.0,
            P10: 0.0,
            P11: INITIAL_DCM_VARIANCE,
            P12: 0.0,
            P13: 0.0,
            P14: 0.0,
            P15: 0.0,
            P20: 0.0,
            P21: 0.0,
            P22: INITIAL_DCM_VARIANCE,
            P23: 0.0,
            P24: 0.0,
            P25: 0.0,
            P30: 0.0,
            P31: 0.0,
            P32: 0.0,
            P33: INITIAL_BIAS_VARIANCE,
            P34: 0.0,
            P35: 0.0,
            P40: 0.0,
            P41: 0.0,
            P42: 0.0,
            P43: 0.0,
            P44: INITIAL_BIAS_VARIANCE,
            P45: 0.0,
            P50: 0.0,
            P51: 0.0,
            P52: 0.0,
            P53: 0.0,
            P54: 0.0,
            P55: INITIAL_BIAS_VARIANCE,
        }
    }

    /// Updates DCMIMU states with gyro (x, y, z), accel (x, y, z),
    /// and dt (seconds) and returns current estimations ({roll; yaw; pitch}).
    pub fn update(
        &mut self,
        gyro: (f32, f32, f32),
        accel: (f32, f32, f32),
        dt: f32,
    ) -> TaitBryanAngles {
        let gx = gyro.0;
        let gy = gyro.1;
        let gz = gyro.2;
        let ax = accel.0;
        let ay = accel.1;
        let az = accel.2;
        // save last state for rotation estimation
        let x_last: [f32; 3] = [self.x0, self.x1, self.x2];
        // state prediction
        let x_0 =
            self.x0 - dt * (gy * self.x2 - gz * self.x1 + self.x1 * self.x5 - self.x2 * self.x4);
        let x_1 =
            self.x1 + dt * (gx * self.x2 - gz * self.x0 + self.x0 * self.x5 - self.x2 * self.x3);
        let x_2 =
            self.x2 - dt * (gx * self.x1 - gy * self.x0 + self.x0 * self.x4 - self.x1 * self.x3);
        let x_3 = self.x3;
        let x_4 = self.x4;
        let x_5 = self.x5;
        // covariance prediction
        let dt2 = dt * dt;
        let P_00 = self.P00
            - dt * (self.P05 * self.x1 - self.P04 * self.x2 - self.P40 * self.x2
                + self.P50 * self.x1
                + self.P02 * (gy - self.x4)
                + self.P20 * (gy - self.x4)
                - self.P01 * (gz - self.x5)
                - self.P10 * (gz - self.x5))
            + dt2
                * (self.q_dcm2
                    - self.x1
                        * (self.P45 * self.x2 - self.P55 * self.x1 - self.P25 * (gy - self.x4)
                            + self.P15 * (gz - self.x5))
                    + self.x2
                        * (self.P44 * self.x2 - self.P54 * self.x1 - self.P24 * (gy - self.x4)
                            + self.P14 * (gz - self.x5))
                    - (gy - self.x4)
                        * (self.P42 * self.x2 - self.P52 * self.x1 - self.P22 * (gy - self.x4)
                            + self.P12 * (gz - self.x5))
                    + (gz - self.x5)
                        * (self.P41 * self.x2 - self.P51 * self.x1 - self.P21 * (gy - self.x4)
                            + self.P11 * (gz - self.x5)));
        let P_01 = self.P01
            + dt * (self.P05 * self.x0 - self.P03 * self.x2 + self.P41 * self.x2
                - self.P51 * self.x1
                + self.P02 * (gx - self.x3)
                - self.P00 * (gz - self.x5)
                - self.P21 * (gy - self.x4)
                + self.P11 * (gz - self.x5))
            + dt2
                * (self.x0
                    * (self.P45 * self.x2 - self.P55 * self.x1 - self.P25 * (gy - self.x4)
                        + self.P15 * (gz - self.x5))
                    - self.x2
                        * (self.P43 * self.x2 - self.P53 * self.x1 - self.P23 * (gy - self.x4)
                            + self.P13 * (gz - self.x5))
                    + (gx - self.x3)
                        * (self.P42 * self.x2 - self.P52 * self.x1 - self.P22 * (gy - self.x4)
                            + self.P12 * (gz - self.x5))
                    - (gz - self.x5)
                        * (self.P40 * self.x2 - self.P50 * self.x1 - self.P20 * (gy - self.x4)
                            + self.P10 * (gz - self.x5)));
        let P_02 = self.P02
            - dt * (self.P04 * self.x0 - self.P03 * self.x1 - self.P42 * self.x2
                + self.P52 * self.x1
                + self.P01 * (gx - self.x3)
                - self.P00 * (gy - self.x4)
                + self.P22 * (gy - self.x4)
                - self.P12 * (gz - self.x5))
            - dt2
                * (self.x0
                    * (self.P44 * self.x2 - self.P54 * self.x1 - self.P24 * (gy - self.x4)
                        + self.P14 * (gz - self.x5))
                    - self.x1
                        * (self.P43 * self.x2 - self.P53 * self.x1 - self.P23 * (gy - self.x4)
                            + self.P13 * (gz - self.x5))
                    + (gx - self.x3)
                        * (self.P41 * self.x2 - self.P51 * self.x1 - self.P21 * (gy - self.x4)
                            + self.P11 * (gz - self.x5))
                    - (gy - self.x4)
                        * (self.P40 * self.x2 - self.P50 * self.x1 - self.P20 * (gy - self.x4)
                            + self.P10 * (gz - self.x5)));
        let P_03 = self.P03
            + dt * (self.P43 * self.x2 - self.P53 * self.x1 - self.P23 * (gy - self.x4)
                + self.P13 * (gz - self.x5));
        let P_04 = self.P04
            + dt * (self.P44 * self.x2 - self.P54 * self.x1 - self.P24 * (gy - self.x4)
                + self.P14 * (gz - self.x5));
        let P_05 = self.P05
            + dt * (self.P45 * self.x2 - self.P55 * self.x1 - self.P25 * (gy - self.x4)
                + self.P15 * (gz - self.x5));
        let P_10 = self.P10
            - dt * (self.P15 * self.x1 - self.P14 * self.x2 + self.P30 * self.x2
                - self.P50 * self.x0
                - self.P20 * (gx - self.x3)
                + self.P12 * (gy - self.x4)
                + self.P00 * (gz - self.x5)
                - self.P11 * (gz - self.x5))
            + dt2
                * (self.x1
                    * (self.P35 * self.x2 - self.P55 * self.x0 - self.P25 * (gx - self.x3)
                        + self.P05 * (gz - self.x5))
                    - self.x2
                        * (self.P34 * self.x2 - self.P54 * self.x0 - self.P24 * (gx - self.x3)
                            + self.P04 * (gz - self.x5))
                    + (gy - self.x4)
                        * (self.P32 * self.x2 - self.P52 * self.x0 - self.P22 * (gx - self.x3)
                            + self.P02 * (gz - self.x5))
                    - (gz - self.x5)
                        * (self.P31 * self.x2 - self.P51 * self.x0 - self.P21 * (gx - self.x3)
                            + self.P01 * (gz - self.x5)));
        let P_11 = self.P11
            + dt * (self.P15 * self.x0 - self.P13 * self.x2 - self.P31 * self.x2
                + self.P51 * self.x0
                + self.P12 * (gx - self.x3)
                + self.P21 * (gx - self.x3)
                - self.P01 * (gz - self.x5)
                - self.P10 * (gz - self.x5))
            + dt2
                * (self.q_dcm2
                    - self.x0
                        * (self.P35 * self.x2 - self.P55 * self.x0 - self.P25 * (gx - self.x3)
                            + self.P05 * (gz - self.x5))
                    + self.x2
                        * (self.P33 * self.x2 - self.P53 * self.x0 - self.P23 * (gx - self.x3)
                            + self.P03 * (gz - self.x5))
                    - (gx - self.x3)
                        * (self.P32 * self.x2 - self.P52 * self.x0 - self.P22 * (gx - self.x3)
                            + self.P02 * (gz - self.x5))
                    + (gz - self.x5)
                        * (self.P30 * self.x2 - self.P50 * self.x0 - self.P20 * (gx - self.x3)
                            + self.P00 * (gz - self.x5)));
        let P_12 = self.P12
            - dt * (self.P14 * self.x0 - self.P13 * self.x1 + self.P32 * self.x2
                - self.P52 * self.x0
                + self.P11 * (gx - self.x3)
                - self.P22 * (gx - self.x3)
                - self.P10 * (gy - self.x4)
                + self.P02 * (gz - self.x5))
            + dt2
                * (self.x0
                    * (self.P34 * self.x2 - self.P54 * self.x0 - self.P24 * (gx - self.x3)
                        + self.P04 * (gz - self.x5))
                    - self.x1
                        * (self.P33 * self.x2 - self.P53 * self.x0 - self.P23 * (gx - self.x3)
                            + self.P03 * (gz - self.x5))
                    + (gx - self.x3)
                        * (self.P31 * self.x2 - self.P51 * self.x0 - self.P21 * (gx - self.x3)
                            + self.P01 * (gz - self.x5))
                    - (gy - self.x4)
                        * (self.P30 * self.x2 - self.P50 * self.x0 - self.P20 * (gx - self.x3)
                            + self.P00 * (gz - self.x5)));
        let P_13 = self.P13
            - dt * (self.P33 * self.x2 - self.P53 * self.x0 - self.P23 * (gx - self.x3)
                + self.P03 * (gz - self.x5));
        let P_14 = self.P14
            - dt * (self.P34 * self.x2 - self.P54 * self.x0 - self.P24 * (gx - self.x3)
                + self.P04 * (gz - self.x5));
        let P_15 = self.P15
            - dt * (self.P35 * self.x2 - self.P55 * self.x0 - self.P25 * (gx - self.x3)
                + self.P05 * (gz - self.x5));
        let P_20 = self.P20
            - dt * (self.P25 * self.x1 - self.P30 * self.x1 + self.P40 * self.x0
                - self.P24 * self.x2
                + self.P10 * (gx - self.x3)
                - self.P00 * (gy - self.x4)
                + self.P22 * (gy - self.x4)
                - self.P21 * (gz - self.x5))
            - dt2
                * (self.x1
                    * (self.P35 * self.x1 - self.P45 * self.x0 - self.P15 * (gx - self.x3)
                        + self.P05 * (gy - self.x4))
                    - self.x2
                        * (self.P34 * self.x1 - self.P44 * self.x0 - self.P14 * (gx - self.x3)
                            + self.P04 * (gy - self.x4))
                    + (gy - self.x4)
                        * (self.P32 * self.x1 - self.P42 * self.x0 - self.P12 * (gx - self.x3)
                            + self.P02 * (gy - self.x4))
                    - (gz - self.x5)
                        * (self.P31 * self.x1 - self.P41 * self.x0 - self.P11 * (gx - self.x3)
                            + self.P01 * (gy - self.x4)));
        let P_21 = self.P21
            + dt * (self.P25 * self.x0 + self.P31 * self.x1
                - self.P41 * self.x0
                - self.P23 * self.x2
                - self.P11 * (gx - self.x3)
                + self.P01 * (gy - self.x4)
                + self.P22 * (gx - self.x3)
                - self.P20 * (gz - self.x5))
            + dt2
                * (self.x0
                    * (self.P35 * self.x1 - self.P45 * self.x0 - self.P15 * (gx - self.x3)
                        + self.P05 * (gy - self.x4))
                    - self.x2
                        * (self.P33 * self.x1 - self.P43 * self.x0 - self.P13 * (gx - self.x3)
                            + self.P03 * (gy - self.x4))
                    + (gx - self.x3)
                        * (self.P32 * self.x1 - self.P42 * self.x0 - self.P12 * (gx - self.x3)
                            + self.P02 * (gy - self.x4))
                    - (gz - self.x5)
                        * (self.P30 * self.x1 - self.P40 * self.x0 - self.P10 * (gx - self.x3)
                            + self.P00 * (gy - self.x4)));
        let P_22 = self.P22
            - dt * (self.P24 * self.x0 - self.P23 * self.x1 - self.P32 * self.x1
                + self.P42 * self.x0
                + self.P12 * (gx - self.x3)
                + self.P21 * (gx - self.x3)
                - self.P02 * (gy - self.x4)
                - self.P20 * (gy - self.x4))
            + dt2
                * (self.q_dcm2
                    - self.x0
                        * (self.P34 * self.x1 - self.P44 * self.x0 - self.P14 * (gx - self.x3)
                            + self.P04 * (gy - self.x4))
                    + self.x1
                        * (self.P33 * self.x1 - self.P43 * self.x0 - self.P13 * (gx - self.x3)
                            + self.P03 * (gy - self.x4))
                    - (gx - self.x3)
                        * (self.P31 * self.x1 - self.P41 * self.x0 - self.P11 * (gx - self.x3)
                            + self.P01 * (gy - self.x4))
                    + (gy - self.x4)
                        * (self.P30 * self.x1 - self.P40 * self.x0 - self.P10 * (gx - self.x3)
                            + self.P00 * (gy - self.x4)));
        let P_23 = self.P23
            + dt * (self.P33 * self.x1 - self.P43 * self.x0 - self.P13 * (gx - self.x3)
                + self.P03 * (gy - self.x4));
        let P_24 = self.P24
            + dt * (self.P34 * self.x1 - self.P44 * self.x0 - self.P14 * (gx - self.x3)
                + self.P04 * (gy - self.x4));
        let P_25 = self.P25
            + dt * (self.P35 * self.x1 - self.P45 * self.x0 - self.P15 * (gx - self.x3)
                + self.P05 * (gy - self.x4));
        let P_30 = self.P30
            - dt * (self.P35 * self.x1 - self.P34 * self.x2 + self.P32 * (gy - self.x4)
                - self.P31 * (gz - self.x5));
        let P_31 = self.P31
            + dt * (self.P35 * self.x0 - self.P33 * self.x2 + self.P32 * (gx - self.x3)
                - self.P30 * (gz - self.x5));
        let P_32 = self.P32
            - dt * (self.P34 * self.x0 - self.P33 * self.x1 + self.P31 * (gx - self.x3)
                - self.P30 * (gy - self.x4));
        let P_33 = self.P33 + dt2 * self.q_gyro_bias2;
        let P_34 = self.P34;
        let P_35 = self.P35;
        let P_40 = self.P40
            - dt * (self.P45 * self.x1 - self.P44 * self.x2 + self.P42 * (gy - self.x4)
                - self.P41 * (gz - self.x5));
        let P_41 = self.P41
            + dt * (self.P45 * self.x0 - self.P43 * self.x2 + self.P42 * (gx - self.x3)
                - self.P40 * (gz - self.x5));
        let P_42 = self.P42
            - dt * (self.P44 * self.x0 - self.P43 * self.x1 + self.P41 * (gx - self.x3)
                - self.P40 * (gy - self.x4));
        let P_43 = self.P43;
        let P_44 = self.P44 + dt2 * self.q_gyro_bias2;
        let P_45 = self.P45;
        let P_50 = self.P50
            - dt * (self.P55 * self.x1 - self.P54 * self.x2 + self.P52 * (gy - self.x4)
                - self.P51 * (gz - self.x5));
        let P_51 = self.P51
            + dt * (self.P55 * self.x0 - self.P53 * self.x2 + self.P52 * (gx - self.x3)
                - self.P50 * (gz - self.x5));
        let P_52 = self.P52
            - dt * (self.P54 * self.x0 - self.P53 * self.x1 + self.P51 * (gx - self.x3)
                - self.P50 * (gy - self.x4));
        let P_53 = self.P53;
        let P_54 = self.P54;
        let P_55 = self.P55 + dt2 * self.q_gyro_bias2;

        // Kalman innovation
        let y0 = ax - self.g0 * x_0;
        let y1 = ay - self.g0 * x_1;
        let y2 = az - self.g0 * x_2;
        let a_len = sqrtf(y0 * y0 + y1 * y1 + y2 * y2);

        let S00 = self.r_acc2 + a_len * self.r_a2 + P_00 * self.g0_2;
        let S01 = P_01 * self.g0_2;
        let S02 = P_02 * self.g0_2;
        let S10 = P_10 * self.g0_2;
        let S11 = self.r_acc2 + a_len * self.r_a2 + P_11 * self.g0_2;
        let S12 = P_12 * self.g0_2;
        let S20 = P_20 * self.g0_2;
        let S21 = P_21 * self.g0_2;
        let S22 = self.r_acc2 + a_len * self.r_a2 + P_22 * self.g0_2;

        // Kalman gain
        let invPart = 1.0
            / (S00 * S11 * S22 - S00 * S12 * S21 - S01 * S10 * S22
                + S01 * S12 * S20
                + S02 * S10 * S21
                - S02 * S11 * S20);
        let K00 = (self.g0
            * (P_02 * S10 * S21 - P_02 * S11 * S20 - P_01 * S10 * S22
                + P_01 * S12 * S20
                + P_00 * S11 * S22
                - P_00 * S12 * S21))
            * invPart;
        let K01 = -(self.g0
            * (P_02 * S00 * S21 - P_02 * S01 * S20 - P_01 * S00 * S22
                + P_01 * S02 * S20
                + P_00 * S01 * S22
                - P_00 * S02 * S21))
            * invPart;
        let K02 = (self.g0
            * (P_02 * S00 * S11 - P_02 * S01 * S10 - P_01 * S00 * S12
                + P_01 * S02 * S10
                + P_00 * S01 * S12
                - P_00 * S02 * S11))
            * invPart;
        let K10 = (self.g0
            * (P_12 * S10 * S21 - P_12 * S11 * S20 - P_11 * S10 * S22
                + P_11 * S12 * S20
                + P_10 * S11 * S22
                - P_10 * S12 * S21))
            * invPart;
        let K11 = -(self.g0
            * (P_12 * S00 * S21 - P_12 * S01 * S20 - P_11 * S00 * S22
                + P_11 * S02 * S20
                + P_10 * S01 * S22
                - P_10 * S02 * S21))
            * invPart;
        let K12 = (self.g0
            * (P_12 * S00 * S11 - P_12 * S01 * S10 - P_11 * S00 * S12
                + P_11 * S02 * S10
                + P_10 * S01 * S12
                - P_10 * S02 * S11))
            * invPart;
        let K20 = (self.g0
            * (P_22 * S10 * S21 - P_22 * S11 * S20 - P_21 * S10 * S22
                + P_21 * S12 * S20
                + P_20 * S11 * S22
                - P_20 * S12 * S21))
            * invPart;
        let K21 = -(self.g0
            * (P_22 * S00 * S21 - P_22 * S01 * S20 - P_21 * S00 * S22
                + P_21 * S02 * S20
                + P_20 * S01 * S22
                - P_20 * S02 * S21))
            * invPart;
        let K22 = (self.g0
            * (P_22 * S00 * S11 - P_22 * S01 * S10 - P_21 * S00 * S12
                + P_21 * S02 * S10
                + P_20 * S01 * S12
                - P_20 * S02 * S11))
            * invPart;
        let K30 = (self.g0
            * (P_32 * S10 * S21 - P_32 * S11 * S20 - P_31 * S10 * S22
                + P_31 * S12 * S20
                + P_30 * S11 * S22
                - P_30 * S12 * S21))
            * invPart;
        let K31 = -(self.g0
            * (P_32 * S00 * S21 - P_32 * S01 * S20 - P_31 * S00 * S22
                + P_31 * S02 * S20
                + P_30 * S01 * S22
                - P_30 * S02 * S21))
            * invPart;
        let K32 = (self.g0
            * (P_32 * S00 * S11 - P_32 * S01 * S10 - P_31 * S00 * S12
                + P_31 * S02 * S10
                + P_30 * S01 * S12
                - P_30 * S02 * S11))
            * invPart;
        let K40 = (self.g0
            * (P_42 * S10 * S21 - P_42 * S11 * S20 - P_41 * S10 * S22
                + P_41 * S12 * S20
                + P_40 * S11 * S22
                - P_40 * S12 * S21))
            * invPart;
        let K41 = -(self.g0
            * (P_42 * S00 * S21 - P_42 * S01 * S20 - P_41 * S00 * S22
                + P_41 * S02 * S20
                + P_40 * S01 * S22
                - P_40 * S02 * S21))
            * invPart;
        let K42 = (self.g0
            * (P_42 * S00 * S11 - P_42 * S01 * S10 - P_41 * S00 * S12
                + P_41 * S02 * S10
                + P_40 * S01 * S12
                - P_40 * S02 * S11))
            * invPart;
        let K50 = (self.g0
            * (P_52 * S10 * S21 - P_52 * S11 * S20 - P_51 * S10 * S22
                + P_51 * S12 * S20
                + P_50 * S11 * S22
                - P_50 * S12 * S21))
            * invPart;
        let K51 = -(self.g0
            * (P_52 * S00 * S21 - P_52 * S01 * S20 - P_51 * S00 * S22
                + P_51 * S02 * S20
                + P_50 * S01 * S22
                - P_50 * S02 * S21))
            * invPart;
        let K52 = (self.g0
            * (P_52 * S00 * S11 - P_52 * S01 * S10 - P_51 * S00 * S12
                + P_51 * S02 * S10
                + P_50 * S01 * S12
                - P_50 * S02 * S11))
            * invPart;

        // update a posteriori
        self.x0 = x_0 + K00 * y0 + K01 * y1 + K02 * y2;
        self.x1 = x_1 + K10 * y0 + K11 * y1 + K12 * y2;
        self.x2 = x_2 + K20 * y0 + K21 * y1 + K22 * y2;
        self.x3 = x_3 + K30 * y0 + K31 * y1 + K32 * y2;
        self.x4 = x_4 + K40 * y0 + K41 * y1 + K42 * y2;
        self.x5 = x_5 + K50 * y0 + K51 * y1 + K52 * y2;

        // update a posteriori covariance
        let r_adab = self.r_acc2 + a_len * self.r_a2;
        let P__00 = P_00
            - self.g0 * (K00 * P_00 * 2.0 + K01 * P_01 + K01 * P_10 + K02 * P_02 + K02 * P_20)
            + (K00 * K00) * r_adab
            + (K01 * K01) * r_adab
            + (K02 * K02) * r_adab
            + self.g0_2
                * (K00 * (K00 * P_00 + K01 * P_10 + K02 * P_20)
                    + K01 * (K00 * P_01 + K01 * P_11 + K02 * P_21)
                    + K02 * (K00 * P_02 + K01 * P_12 + K02 * P_22));
        let P__01 = P_01
            - self.g0
                * (K00 * P_01 + K01 * P_11 + K02 * P_21 + K10 * P_00 + K11 * P_01 + K12 * P_02)
            + self.g0_2
                * (K10 * (K00 * P_00 + K01 * P_10 + K02 * P_20)
                    + K11 * (K00 * P_01 + K01 * P_11 + K02 * P_21)
                    + K12 * (K00 * P_02 + K01 * P_12 + K02 * P_22))
            + K00 * K10 * r_adab
            + K01 * K11 * r_adab
            + K02 * K12 * r_adab;
        let P__02 = P_02
            - self.g0
                * (K00 * P_02 + K01 * P_12 + K02 * P_22 + K20 * P_00 + K21 * P_01 + K22 * P_02)
            + self.g0_2
                * (K20 * (K00 * P_00 + K01 * P_10 + K02 * P_20)
                    + K21 * (K00 * P_01 + K01 * P_11 + K02 * P_21)
                    + K22 * (K00 * P_02 + K01 * P_12 + K02 * P_22))
            + K00 * K20 * r_adab
            + K01 * K21 * r_adab
            + K02 * K22 * r_adab;
        let P__03 = P_03
            - self.g0
                * (K00 * P_03 + K01 * P_13 + K02 * P_23 + K30 * P_00 + K31 * P_01 + K32 * P_02)
            + self.g0_2
                * (K30 * (K00 * P_00 + K01 * P_10 + K02 * P_20)
                    + K31 * (K00 * P_01 + K01 * P_11 + K02 * P_21)
                    + K32 * (K00 * P_02 + K01 * P_12 + K02 * P_22))
            + K00 * K30 * r_adab
            + K01 * K31 * r_adab
            + K02 * K32 * r_adab;
        let P__04 = P_04
            - self.g0
                * (K00 * P_04 + K01 * P_14 + K02 * P_24 + K40 * P_00 + K41 * P_01 + K42 * P_02)
            + self.g0_2
                * (K40 * (K00 * P_00 + K01 * P_10 + K02 * P_20)
                    + K41 * (K00 * P_01 + K01 * P_11 + K02 * P_21)
                    + K42 * (K00 * P_02 + K01 * P_12 + K02 * P_22))
            + K00 * K40 * r_adab
            + K01 * K41 * r_adab
            + K02 * K42 * r_adab;
        let P__05 = P_05
            - self.g0
                * (K00 * P_05 + K01 * P_15 + K02 * P_25 + K50 * P_00 + K51 * P_01 + K52 * P_02)
            + self.g0_2
                * (K50 * (K00 * P_00 + K01 * P_10 + K02 * P_20)
                    + K51 * (K00 * P_01 + K01 * P_11 + K02 * P_21)
                    + K52 * (K00 * P_02 + K01 * P_12 + K02 * P_22))
            + K00 * K50 * r_adab
            + K01 * K51 * r_adab
            + K02 * K52 * r_adab;
        let P__10 = P_10
            - self.g0
                * (K00 * P_10 + K01 * P_11 + K02 * P_12 + K10 * P_00 + K11 * P_10 + K12 * P_20)
            + self.g0_2
                * (K00 * (K10 * P_00 + K11 * P_10 + K12 * P_20)
                    + K01 * (K10 * P_01 + K11 * P_11 + K12 * P_21)
                    + K02 * (K10 * P_02 + K11 * P_12 + K12 * P_22))
            + K00 * K10 * r_adab
            + K01 * K11 * r_adab
            + K02 * K12 * r_adab;
        let P__11 = P_11
            - self.g0 * (K10 * P_01 + K10 * P_10 + K11 * P_11 * 2.0 + K12 * P_12 + K12 * P_21)
            + (K10 * K10) * r_adab
            + (K11 * K11) * r_adab
            + (K12 * K12) * r_adab
            + self.g0_2
                * (K10 * (K10 * P_00 + K11 * P_10 + K12 * P_20)
                    + K11 * (K10 * P_01 + K11 * P_11 + K12 * P_21)
                    + K12 * (K10 * P_02 + K11 * P_12 + K12 * P_22));
        let P__12 = P_12
            - self.g0
                * (K10 * P_02 + K11 * P_12 + K12 * P_22 + K20 * P_10 + K21 * P_11 + K22 * P_12)
            + self.g0_2
                * (K20 * (K10 * P_00 + K11 * P_10 + K12 * P_20)
                    + K21 * (K10 * P_01 + K11 * P_11 + K12 * P_21)
                    + K22 * (K10 * P_02 + K11 * P_12 + K12 * P_22))
            + K10 * K20 * r_adab
            + K11 * K21 * r_adab
            + K12 * K22 * r_adab;
        let P__13 = P_13
            - self.g0
                * (K10 * P_03 + K11 * P_13 + K12 * P_23 + K30 * P_10 + K31 * P_11 + K32 * P_12)
            + self.g0_2
                * (K30 * (K10 * P_00 + K11 * P_10 + K12 * P_20)
                    + K31 * (K10 * P_01 + K11 * P_11 + K12 * P_21)
                    + K32 * (K10 * P_02 + K11 * P_12 + K12 * P_22))
            + K10 * K30 * r_adab
            + K11 * K31 * r_adab
            + K12 * K32 * r_adab;
        let P__14 = P_14
            - self.g0
                * (K10 * P_04 + K11 * P_14 + K12 * P_24 + K40 * P_10 + K41 * P_11 + K42 * P_12)
            + self.g0_2
                * (K40 * (K10 * P_00 + K11 * P_10 + K12 * P_20)
                    + K41 * (K10 * P_01 + K11 * P_11 + K12 * P_21)
                    + K42 * (K10 * P_02 + K11 * P_12 + K12 * P_22))
            + K10 * K40 * r_adab
            + K11 * K41 * r_adab
            + K12 * K42 * r_adab;
        let P__15 = P_15
            - self.g0
                * (K10 * P_05 + K11 * P_15 + K12 * P_25 + K50 * P_10 + K51 * P_11 + K52 * P_12)
            + self.g0_2
                * (K50 * (K10 * P_00 + K11 * P_10 + K12 * P_20)
                    + K51 * (K10 * P_01 + K11 * P_11 + K12 * P_21)
                    + K52 * (K10 * P_02 + K11 * P_12 + K12 * P_22))
            + K10 * K50 * r_adab
            + K11 * K51 * r_adab
            + K12 * K52 * r_adab;
        let P__20 = P_20
            - self.g0
                * (K00 * P_20 + K01 * P_21 + K02 * P_22 + K20 * P_00 + K21 * P_10 + K22 * P_20)
            + self.g0_2
                * (K00 * (K20 * P_00 + K21 * P_10 + K22 * P_20)
                    + K01 * (K20 * P_01 + K21 * P_11 + K22 * P_21)
                    + K02 * (K20 * P_02 + K21 * P_12 + K22 * P_22))
            + K00 * K20 * r_adab
            + K01 * K21 * r_adab
            + K02 * K22 * r_adab;
        let P__21 = P_21
            - self.g0
                * (K10 * P_20 + K11 * P_21 + K12 * P_22 + K20 * P_01 + K21 * P_11 + K22 * P_21)
            + self.g0_2
                * (K10 * (K20 * P_00 + K21 * P_10 + K22 * P_20)
                    + K11 * (K20 * P_01 + K21 * P_11 + K22 * P_21)
                    + K12 * (K20 * P_02 + K21 * P_12 + K22 * P_22))
            + K10 * K20 * r_adab
            + K11 * K21 * r_adab
            + K12 * K22 * r_adab;
        let P__22 = P_22
            - self.g0 * (K20 * P_02 + K20 * P_20 + K21 * P_12 + K21 * P_21 + K22 * P_22 * 2.0)
            + (K20 * K20) * r_adab
            + (K21 * K21) * r_adab
            + (K22 * K22) * r_adab
            + self.g0_2
                * (K20 * (K20 * P_00 + K21 * P_10 + K22 * P_20)
                    + K21 * (K20 * P_01 + K21 * P_11 + K22 * P_21)
                    + K22 * (K20 * P_02 + K21 * P_12 + K22 * P_22));
        let P__23 = P_23
            - self.g0
                * (K20 * P_03 + K21 * P_13 + K22 * P_23 + K30 * P_20 + K31 * P_21 + K32 * P_22)
            + self.g0_2
                * (K30 * (K20 * P_00 + K21 * P_10 + K22 * P_20)
                    + K31 * (K20 * P_01 + K21 * P_11 + K22 * P_21)
                    + K32 * (K20 * P_02 + K21 * P_12 + K22 * P_22))
            + K20 * K30 * r_adab
            + K21 * K31 * r_adab
            + K22 * K32 * r_adab;
        let P__24 = P_24
            - self.g0
                * (K20 * P_04 + K21 * P_14 + K22 * P_24 + K40 * P_20 + K41 * P_21 + K42 * P_22)
            + self.g0_2
                * (K40 * (K20 * P_00 + K21 * P_10 + K22 * P_20)
                    + K41 * (K20 * P_01 + K21 * P_11 + K22 * P_21)
                    + K42 * (K20 * P_02 + K21 * P_12 + K22 * P_22))
            + K20 * K40 * r_adab
            + K21 * K41 * r_adab
            + K22 * K42 * r_adab;
        let P__25 = P_25
            - self.g0
                * (K20 * P_05 + K21 * P_15 + K22 * P_25 + K50 * P_20 + K51 * P_21 + K52 * P_22)
            + self.g0_2
                * (K50 * (K20 * P_00 + K21 * P_10 + K22 * P_20)
                    + K51 * (K20 * P_01 + K21 * P_11 + K22 * P_21)
                    + K52 * (K20 * P_02 + K21 * P_12 + K22 * P_22))
            + K20 * K50 * r_adab
            + K21 * K51 * r_adab
            + K22 * K52 * r_adab;
        let P__30 = P_30
            - self.g0
                * (K00 * P_30 + K01 * P_31 + K02 * P_32 + K30 * P_00 + K31 * P_10 + K32 * P_20)
            + self.g0_2
                * (K00 * (K30 * P_00 + K31 * P_10 + K32 * P_20)
                    + K01 * (K30 * P_01 + K31 * P_11 + K32 * P_21)
                    + K02 * (K30 * P_02 + K31 * P_12 + K32 * P_22))
            + K00 * K30 * r_adab
            + K01 * K31 * r_adab
            + K02 * K32 * r_adab;
        let P__31 = P_31
            - self.g0
                * (K10 * P_30 + K11 * P_31 + K12 * P_32 + K30 * P_01 + K31 * P_11 + K32 * P_21)
            + self.g0_2
                * (K10 * (K30 * P_00 + K31 * P_10 + K32 * P_20)
                    + K11 * (K30 * P_01 + K31 * P_11 + K32 * P_21)
                    + K12 * (K30 * P_02 + K31 * P_12 + K32 * P_22))
            + K10 * K30 * r_adab
            + K11 * K31 * r_adab
            + K12 * K32 * r_adab;
        let P__32 = P_32
            - self.g0
                * (K20 * P_30 + K21 * P_31 + K22 * P_32 + K30 * P_02 + K31 * P_12 + K32 * P_22)
            + self.g0_2
                * (K20 * (K30 * P_00 + K31 * P_10 + K32 * P_20)
                    + K21 * (K30 * P_01 + K31 * P_11 + K32 * P_21)
                    + K22 * (K30 * P_02 + K31 * P_12 + K32 * P_22))
            + K20 * K30 * r_adab
            + K21 * K31 * r_adab
            + K22 * K32 * r_adab;
        let P__33 = P_33
            - self.g0
                * (K30 * P_03 + K31 * P_13 + K30 * P_30 + K31 * P_31 + K32 * P_23 + K32 * P_32)
            + (K30 * K30) * r_adab
            + (K31 * K31) * r_adab
            + (K32 * K32) * r_adab
            + self.g0_2
                * (K30 * (K30 * P_00 + K31 * P_10 + K32 * P_20)
                    + K31 * (K30 * P_01 + K31 * P_11 + K32 * P_21)
                    + K32 * (K30 * P_02 + K31 * P_12 + K32 * P_22));
        let P__34 = P_34
            - self.g0
                * (K30 * P_04 + K31 * P_14 + K32 * P_24 + K40 * P_30 + K41 * P_31 + K42 * P_32)
            + self.g0_2
                * (K40 * (K30 * P_00 + K31 * P_10 + K32 * P_20)
                    + K41 * (K30 * P_01 + K31 * P_11 + K32 * P_21)
                    + K42 * (K30 * P_02 + K31 * P_12 + K32 * P_22))
            + K30 * K40 * r_adab
            + K31 * K41 * r_adab
            + K32 * K42 * r_adab;
        let P__35 = P_35
            - self.g0
                * (K30 * P_05 + K31 * P_15 + K32 * P_25 + K50 * P_30 + K51 * P_31 + K52 * P_32)
            + self.g0_2
                * (K50 * (K30 * P_00 + K31 * P_10 + K32 * P_20)
                    + K51 * (K30 * P_01 + K31 * P_11 + K32 * P_21)
                    + K52 * (K30 * P_02 + K31 * P_12 + K32 * P_22))
            + K30 * K50 * r_adab
            + K31 * K51 * r_adab
            + K32 * K52 * r_adab;
        let P__40 = P_40
            - self.g0
                * (K00 * P_40 + K01 * P_41 + K02 * P_42 + K40 * P_00 + K41 * P_10 + K42 * P_20)
            + self.g0_2
                * (K00 * (K40 * P_00 + K41 * P_10 + K42 * P_20)
                    + K01 * (K40 * P_01 + K41 * P_11 + K42 * P_21)
                    + K02 * (K40 * P_02 + K41 * P_12 + K42 * P_22))
            + K00 * K40 * r_adab
            + K01 * K41 * r_adab
            + K02 * K42 * r_adab;
        let P__41 = P_41
            - self.g0
                * (K10 * P_40 + K11 * P_41 + K12 * P_42 + K40 * P_01 + K41 * P_11 + K42 * P_21)
            + self.g0_2
                * (K10 * (K40 * P_00 + K41 * P_10 + K42 * P_20)
                    + K11 * (K40 * P_01 + K41 * P_11 + K42 * P_21)
                    + K12 * (K40 * P_02 + K41 * P_12 + K42 * P_22))
            + K10 * K40 * r_adab
            + K11 * K41 * r_adab
            + K12 * K42 * r_adab;
        let P__42 = P_42
            - self.g0
                * (K20 * P_40 + K21 * P_41 + K22 * P_42 + K40 * P_02 + K41 * P_12 + K42 * P_22)
            + self.g0_2
                * (K20 * (K40 * P_00 + K41 * P_10 + K42 * P_20)
                    + K21 * (K40 * P_01 + K41 * P_11 + K42 * P_21)
                    + K22 * (K40 * P_02 + K41 * P_12 + K42 * P_22))
            + K20 * K40 * r_adab
            + K21 * K41 * r_adab
            + K22 * K42 * r_adab;
        let P__43 = P_43
            - self.g0
                * (K30 * P_40 + K31 * P_41 + K32 * P_42 + K40 * P_03 + K41 * P_13 + K42 * P_23)
            + self.g0_2
                * (K30 * (K40 * P_00 + K41 * P_10 + K42 * P_20)
                    + K31 * (K40 * P_01 + K41 * P_11 + K42 * P_21)
                    + K32 * (K40 * P_02 + K41 * P_12 + K42 * P_22))
            + K30 * K40 * r_adab
            + K31 * K41 * r_adab
            + K32 * K42 * r_adab;
        let P__44 = P_44
            - self.g0
                * (K40 * P_04 + K41 * P_14 + K40 * P_40 + K42 * P_24 + K41 * P_41 + K42 * P_42)
            + (K40 * K40) * r_adab
            + (K41 * K41) * r_adab
            + (K42 * K42) * r_adab
            + self.g0_2
                * (K40 * (K40 * P_00 + K41 * P_10 + K42 * P_20)
                    + K41 * (K40 * P_01 + K41 * P_11 + K42 * P_21)
                    + K42 * (K40 * P_02 + K41 * P_12 + K42 * P_22));
        let P__45 = P_45
            - self.g0
                * (K40 * P_05 + K41 * P_15 + K42 * P_25 + K50 * P_40 + K51 * P_41 + K52 * P_42)
            + self.g0_2
                * (K50 * (K40 * P_00 + K41 * P_10 + K42 * P_20)
                    + K51 * (K40 * P_01 + K41 * P_11 + K42 * P_21)
                    + K52 * (K40 * P_02 + K41 * P_12 + K42 * P_22))
            + K40 * K50 * r_adab
            + K41 * K51 * r_adab
            + K42 * K52 * r_adab;
        let P__50 = P_50
            - self.g0
                * (K00 * P_50 + K01 * P_51 + K02 * P_52 + K50 * P_00 + K51 * P_10 + K52 * P_20)
            + self.g0_2
                * (K00 * (K50 * P_00 + K51 * P_10 + K52 * P_20)
                    + K01 * (K50 * P_01 + K51 * P_11 + K52 * P_21)
                    + K02 * (K50 * P_02 + K51 * P_12 + K52 * P_22))
            + K00 * K50 * r_adab
            + K01 * K51 * r_adab
            + K02 * K52 * r_adab;
        let P__51 = P_51
            - self.g0
                * (K10 * P_50 + K11 * P_51 + K12 * P_52 + K50 * P_01 + K51 * P_11 + K52 * P_21)
            + self.g0_2
                * (K10 * (K50 * P_00 + K51 * P_10 + K52 * P_20)
                    + K11 * (K50 * P_01 + K51 * P_11 + K52 * P_21)
                    + K12 * (K50 * P_02 + K51 * P_12 + K52 * P_22))
            + K10 * K50 * r_adab
            + K11 * K51 * r_adab
            + K12 * K52 * r_adab;
        let P__52 = P_52
            - self.g0
                * (K20 * P_50 + K21 * P_51 + K22 * P_52 + K50 * P_02 + K51 * P_12 + K52 * P_22)
            + self.g0_2
                * (K20 * (K50 * P_00 + K51 * P_10 + K52 * P_20)
                    + K21 * (K50 * P_01 + K51 * P_11 + K52 * P_21)
                    + K22 * (K50 * P_02 + K51 * P_12 + K52 * P_22))
            + K20 * K50 * r_adab
            + K21 * K51 * r_adab
            + K22 * K52 * r_adab;
        let P__53 = P_53
            - self.g0
                * (K30 * P_50 + K31 * P_51 + K32 * P_52 + K50 * P_03 + K51 * P_13 + K52 * P_23)
            + self.g0_2
                * (K30 * (K50 * P_00 + K51 * P_10 + K52 * P_20)
                    + K31 * (K50 * P_01 + K51 * P_11 + K52 * P_21)
                    + K32 * (K50 * P_02 + K51 * P_12 + K52 * P_22))
            + K30 * K50 * r_adab
            + K31 * K51 * r_adab
            + K32 * K52 * r_adab;
        let P__54 = P_54
            - self.g0
                * (K40 * P_50 + K41 * P_51 + K42 * P_52 + K50 * P_04 + K51 * P_14 + K52 * P_24)
            + self.g0_2
                * (K40 * (K50 * P_00 + K51 * P_10 + K52 * P_20)
                    + K41 * (K50 * P_01 + K51 * P_11 + K52 * P_21)
                    + K42 * (K50 * P_02 + K51 * P_12 + K52 * P_22))
            + K40 * K50 * r_adab
            + K41 * K51 * r_adab
            + K42 * K52 * r_adab;
        let P__55 = P_55
            - self.g0
                * (K50 * P_05 + K51 * P_15 + K52 * P_25 + K50 * P_50 + K51 * P_51 + K52 * P_52)
            + (K50 * K50) * r_adab
            + (K51 * K51) * r_adab
            + (K52 * K52) * r_adab
            + self.g0_2
                * (K50 * (K50 * P_00 + K51 * P_10 + K52 * P_20)
                    + K51 * (K50 * P_01 + K51 * P_11 + K52 * P_21)
                    + K52 * (K50 * P_02 + K51 * P_12 + K52 * P_22));

        let len = sqrtf(self.x0 * self.x0 + self.x1 * self.x1 + self.x2 * self.x2);
        let invlen3 = 1.0 / (len * len * len);
        let invlen32 = invlen3 * invlen3;

        let x1_x2 = self.x1 * self.x1 + self.x2 * self.x2;
        let x0_x2 = self.x0 * self.x0 + self.x2 * self.x2;
        let x0_x1 = self.x0 * self.x0 + self.x1 * self.x1;

        // normalized a posteriori covariance
        self.P00 = invlen32
            * (-x1_x2 * (-P__00 * x1_x2 + P__10 * self.x0 * self.x1 + P__20 * self.x0 * self.x2)
                + self.x0
                    * self.x1
                    * (-P__01 * x1_x2 + P__11 * self.x0 * self.x1 + P__21 * self.x0 * self.x2)
                + self.x0
                    * self.x2
                    * (-P__02 * x1_x2 + P__12 * self.x0 * self.x1 + P__22 * self.x0 * self.x2));
        self.P01 = invlen32
            * (-x0_x2 * (-P__01 * x1_x2 + P__11 * self.x0 * self.x1 + P__21 * self.x0 * self.x2)
                + self.x0
                    * self.x1
                    * (-P__00 * x1_x2 + P__10 * self.x0 * self.x1 + P__20 * self.x0 * self.x2)
                + self.x1
                    * self.x2
                    * (-P__02 * x1_x2 + P__12 * self.x0 * self.x1 + P__22 * self.x0 * self.x2));
        self.P02 = invlen32
            * (-x0_x1 * (-P__02 * x1_x2 + P__12 * self.x0 * self.x1 + P__22 * self.x0 * self.x2)
                + self.x0
                    * self.x2
                    * (-P__00 * x1_x2 + P__10 * self.x0 * self.x1 + P__20 * self.x0 * self.x2)
                + self.x1
                    * self.x2
                    * (-P__01 * x1_x2 + P__11 * self.x0 * self.x1 + P__21 * self.x0 * self.x2));
        self.P03 =
            -invlen3 * (-P__03 * x1_x2 + P__13 * self.x0 * self.x1 + P__23 * self.x0 * self.x2);
        self.P04 =
            -invlen3 * (-P__04 * x1_x2 + P__14 * self.x0 * self.x1 + P__24 * self.x0 * self.x2);
        self.P05 =
            -invlen3 * (-P__05 * x1_x2 + P__15 * self.x0 * self.x1 + P__25 * self.x0 * self.x2);
        self.P10 = invlen32
            * (-x1_x2 * (-P__10 * x0_x2 + P__00 * self.x0 * self.x1 + P__20 * self.x1 * self.x2)
                + self.x0
                    * self.x1
                    * (-P__11 * x0_x2 + P__01 * self.x0 * self.x1 + P__21 * self.x1 * self.x2)
                + self.x0
                    * self.x2
                    * (-P__12 * x0_x2 + P__02 * self.x0 * self.x1 + P__22 * self.x1 * self.x2));
        self.P11 = invlen32
            * (-x0_x2 * (-P__11 * x0_x2 + P__01 * self.x0 * self.x1 + P__21 * self.x1 * self.x2)
                + self.x0
                    * self.x1
                    * (-P__10 * x0_x2 + P__00 * self.x0 * self.x1 + P__20 * self.x1 * self.x2)
                + self.x1
                    * self.x2
                    * (-P__12 * x0_x2 + P__02 * self.x0 * self.x1 + P__22 * self.x1 * self.x2));
        self.P12 = invlen32
            * (-x0_x1 * (-P__12 * x0_x2 + P__02 * self.x0 * self.x1 + P__22 * self.x1 * self.x2)
                + self.x0
                    * self.x2
                    * (-P__10 * x0_x2 + P__00 * self.x0 * self.x1 + P__20 * self.x1 * self.x2)
                + self.x1
                    * self.x2
                    * (-P__11 * x0_x2 + P__01 * self.x0 * self.x1 + P__21 * self.x1 * self.x2));
        self.P13 =
            -invlen3 * (-P__13 * x0_x2 + P__03 * self.x0 * self.x1 + P__23 * self.x1 * self.x2);
        self.P14 =
            -invlen3 * (-P__14 * x0_x2 + P__04 * self.x0 * self.x1 + P__24 * self.x1 * self.x2);
        self.P15 =
            -invlen3 * (-P__15 * x0_x2 + P__05 * self.x0 * self.x1 + P__25 * self.x1 * self.x2);
        self.P20 = invlen32
            * (-x1_x2 * (-P__20 * x0_x1 + P__00 * self.x0 * self.x2 + P__10 * self.x1 * self.x2)
                + self.x0
                    * self.x1
                    * (-P__21 * x0_x1 + P__01 * self.x0 * self.x2 + P__11 * self.x1 * self.x2)
                + self.x0
                    * self.x2
                    * (-P__22 * x0_x1 + P__02 * self.x0 * self.x2 + P__12 * self.x1 * self.x2));
        self.P21 = invlen32
            * (-x0_x2 * (-P__21 * x0_x1 + P__01 * self.x0 * self.x2 + P__11 * self.x1 * self.x2)
                + self.x0
                    * self.x1
                    * (-P__20 * x0_x1 + P__00 * self.x0 * self.x2 + P__10 * self.x1 * self.x2)
                + self.x1
                    * self.x2
                    * (-P__22 * x0_x1 + P__02 * self.x0 * self.x2 + P__12 * self.x1 * self.x2));
        self.P22 = invlen32
            * (-x0_x1 * (-P__22 * x0_x1 + P__02 * self.x0 * self.x2 + P__12 * self.x1 * self.x2)
                + self.x0
                    * self.x2
                    * (-P__20 * x0_x1 + P__00 * self.x0 * self.x2 + P__10 * self.x1 * self.x2)
                + self.x1
                    * self.x2
                    * (-P__21 * x0_x1 + P__01 * self.x0 * self.x2 + P__11 * self.x1 * self.x2));
        self.P23 =
            -invlen3 * (-P__23 * x0_x1 + P__03 * self.x0 * self.x2 + P__13 * self.x1 * self.x2);
        self.P24 =
            -invlen3 * (-P__24 * x0_x1 + P__04 * self.x0 * self.x2 + P__14 * self.x1 * self.x2);
        self.P25 =
            -invlen3 * (-P__25 * x0_x1 + P__05 * self.x0 * self.x2 + P__15 * self.x1 * self.x2);
        self.P30 =
            -invlen3 * (-P__30 * x1_x2 + P__31 * self.x0 * self.x1 + P__32 * self.x0 * self.x2);
        self.P31 =
            -invlen3 * (-P__31 * x0_x2 + P__30 * self.x0 * self.x1 + P__32 * self.x1 * self.x2);
        self.P32 =
            -invlen3 * (-P__32 * x0_x1 + P__30 * self.x0 * self.x2 + P__31 * self.x1 * self.x2);
        self.P33 = P__33;
        self.P34 = P__34;
        self.P35 = P__35;
        self.P40 =
            -invlen3 * (-P__40 * x1_x2 + P__41 * self.x0 * self.x1 + P__42 * self.x0 * self.x2);
        self.P41 =
            -invlen3 * (-P__41 * x0_x2 + P__40 * self.x0 * self.x1 + P__42 * self.x1 * self.x2);
        self.P42 =
            -invlen3 * (-P__42 * x0_x1 + P__40 * self.x0 * self.x2 + P__41 * self.x1 * self.x2);
        self.P43 = P__43;
        self.P44 = P__44;
        self.P45 = P__45;
        self.P50 =
            -invlen3 * (-P__50 * x1_x2 + P__51 * self.x0 * self.x1 + P__52 * self.x0 * self.x2);
        self.P51 =
            -invlen3 * (-P__51 * x0_x2 + P__50 * self.x0 * self.x1 + P__52 * self.x1 * self.x2);
        self.P52 =
            -invlen3 * (-P__52 * x0_x1 + P__50 * self.x0 * self.x2 + P__51 * self.x1 * self.x2);
        self.P53 = P__53;
        self.P54 = P__54;
        self.P55 = P__55;
        // normalized a posteriori state
        self.x0 = self.x0 / len;
        self.x1 = self.x1 / len;
        self.x2 = self.x2 / len;
        // compute Euler angles
        let u_nb1 = gy - self.x4;
        let u_nb2 = gz - self.x5;
        let cy = cosf(self.yaw);
        let sy = sinf(self.yaw);
        let d = sqrtf(x_last[1] * x_last[1] + x_last[2] * x_last[2]);
        let d_inv = 1.0 / d;
        // compute needed parts of rotation matrix R (state and angle based version, equivalent with the commented version above)
        let R11 = cy * d;
        let R12 = -(x_last[2] * sy + x_last[0] * x_last[1] * cy) * d_inv;
        let R13 = (x_last[1] * sy - x_last[0] * x_last[2] * cy) * d_inv;
        let R21 = sy * d;
        let R22 = (x_last[2] * cy - x_last[0] * x_last[1] * sy) * d_inv;
        let R23 = -(x_last[1] * cy + x_last[0] * x_last[2] * sy) * d_inv;

        // update needed parts of R for yaw computation
        let R11_new = R11 + dt * (u_nb2 * R12 - u_nb1 * R13);
        let R21_new = R21 + dt * (u_nb2 * R22 - u_nb1 * R23);

        self.yaw = atan2f(R21_new, R11_new);
        self.pitch = asinf(-self.x0);
        self.roll = atan2f(self.x1, self.x2);

        // save the estimated non-gravitational acceleration
        self.a0 = ax - self.x0 * self.g0;
        self.a1 = ay - self.x1 * self.g0;
        self.a2 = az - self.x2 * self.g0;

        self.all()
    }

    /// Returns all moments (yaw, roll, pitch)
    pub fn all(&self) -> TaitBryanAngles {
        TaitBryanAngles {
            yaw: self.yaw,
            roll: self.roll,
            pitch: self.pitch,
        }
    }

    /// Yaw
    pub fn yaw(&self) -> f32 {
        self.yaw
    }

    /// Pitch
    pub fn pitch(&self) -> f32 {
        self.pitch
    }

    /// Roll
    pub fn roll(&self) -> f32 {
        self.roll
    }
}

/// Represents three dimensions:
///  * yaw, nose left or right about an axis running up and down;
///  * pitch, nose up or down about an axis running from wing to wing;
///  * roll, rotation about an axis running from nose to tail.
/// The axes are alternatively designated as
/// vertical, transverse, and longitudinal respectively.
/// See https://en.wikipedia.org/wiki/Euler_angles#Tait%E2%80%93Bryan_angles
/// and https://en.wikipedia.org/wiki/Aircraft_principal_axes.
#[derive(Debug, Clone, Copy)]
pub struct TaitBryanAngles {
    pub yaw: f32,
    pub pitch: f32,
    pub roll: f32,
}
