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
	P00: f32, P01: f32, P02: f32, P03: f32, P04: f32, P05: f32,
	P10: f32, P11: f32, P12: f32, P13: f32, P14: f32, P15: f32,
	P20: f32, P21: f32, P22: f32, P23: f32, P24: f32, P25: f32,
	P30: f32, P31: f32, P32: f32, P33: f32, P34: f32, P35: f32,
	P40: f32, P41: f32, P42: f32, P43: f32, P44: f32, P45: f32,
	P50: f32, P51: f32, P52: f32, P53: f32, P54: f32, P55: f32,
}

pub const GRAVITY: f32 = 9.81;
const InitialDCMVariance: f32 = 1.0;
const InitialBiasVariance: f32 = 0.1 * 0.1;

impl DCMIMU {
    pub fn new() -> Self {
        DCMIMU {
            g0: GRAVITY,
            g0_2: GRAVITY * GRAVITY,
            x0: 0.0, x1: 0.0, x2: 1.0, x3: 0.0, x4: 0.0, x5: 0.0,
            q_dcm2: 0.1 * 0.1,
            q_gyro_bias2: 0.0001 * 0.0001,
            r_acc2: 0.5 * 0.5,
            r_a2: 10.0 * 10.0,
            a0: 0.0, a1: 0.0, a2: 0.0,
            yaw: 0.0, pitch: 0.0, roll: 0.0,
            P00: InitialDCMVariance, P01:0.0, P02: 0.0, P03: 0.0, P04: 0.0,
            P05: 0.0, P10: 0.0, P11: InitialDCMVariance, P12: 0.0, P13: 0.0,
            P14: 0.0, P15: 0.0, P20: 0.0, P21: 0.0, P22: InitialDCMVariance,
            P23: 0.0, P24: 0.0, P25: 0.0, P30: 0.0, P31: 0.0, P32: 0.0,
            P33: InitialBiasVariance, P34: 0.0, P35: 0.0, P40: 0.0, P41: 0.0,
            P42: 0.0, P43: 0.0, P44: InitialBiasVariance, P45: 0.0, P50: 0.0,
            P51: 0.0, P52: 0.0, P53: 0.0, P54: 0.0, P55: InitialBiasVariance,

        }
    }
}