//! Spectrum containers: the channel/counts vectors plus the gamma header.
//!
//! [`Spectrum`] mirrors the upstream pulse-height container (name, channel
//! numbering, channel labels, energy bins, counts). [`GammaSpectrum`] adds
//! the acquisition header (times, detector id, calibration fits, file name).

use crate::{calib::energy_bins, Error};

/// Pulse-height spectrum container.
#[derive(Debug, Clone, PartialEq)]
pub struct Spectrum {
    /// Spectrum name.
    pub spec_name: String,
    /// First channel number as declared by the file (labels may still be
    /// 0-based positional; see the E-index quirks page).
    pub start_chan_num: i64,
    /// Declared number of channels.
    pub num_channels: usize,
    /// Channel labels (dollar format: positional `0..len`; plain: from file).
    pub channels: Vec<f64>,
    /// Energy per channel (same units as the calibration fit).
    pub ebin: Vec<f64>,
    /// Counts per channel.
    pub counts: Vec<f64>,
}

impl Spectrum {
    /// Empty spectrum with the given name.
    pub fn new(spec_name: &str) -> Self {
        Self {
            spec_name: spec_name.to_string(),
            start_chan_num: 0,
            num_channels: 0,
            channels: Vec::new(),
            ebin: Vec::new(),
            counts: Vec::new(),
        }
    }

    /// Number of channels held.
    pub fn len(&self) -> usize {
        self.counts.len()
    }

    /// Whether any channels are held.
    pub fn is_empty(&self) -> bool {
        self.counts.is_empty()
    }

    /// Largest channel label, for the E3–E5 range checks.
    pub fn max_channel(&self) -> Option<f64> {
        self.channels
            .iter()
            .copied()
            .reduce(|a, b| if a >= b { a } else { b })
    }
}

/// Gamma-ray spectrum: a [`Spectrum`] plus the acquisition header.
#[derive(Debug, Clone, PartialEq)]
pub struct GammaSpectrum {
    /// Channel/counts container.
    pub spectrum: Spectrum,
    /// Real (clock) time.
    pub real_time: f64,
    /// Live time.
    pub live_time: f64,
    /// Detector id.
    pub det_id: String,
    /// Detector description.
    pub det_descp: String,
    /// Acquisition start date (verbatim file text).
    pub start_date: String,
    /// Acquisition start time (verbatim file text).
    pub start_time: String,
    /// Energy calibration fit `[a0, a1, a2]` for E6.
    pub calib_e_fit: Vec<f64>,
    /// FWHM calibration fit (parsed, never evaluated).
    pub calib_fwhm_fit: Vec<f64>,
    /// Source file name (empty when parsed from text).
    pub file_name: String,
}

impl GammaSpectrum {
    /// Empty spectrum with default header values.
    pub fn new() -> Self {
        Self {
            spectrum: Spectrum::new(""),
            real_time: 0.0,
            live_time: 0.0,
            det_id: String::new(),
            det_descp: String::new(),
            start_date: String::new(),
            start_time: String::new(),
            calib_e_fit: Vec::new(),
            calib_fwhm_fit: Vec::new(),
            file_name: String::new(),
        }
    }

    /// Dead time (`real_time - live_time`).
    pub fn dead_time(&self) -> f64 {
        self.real_time - self.live_time
    }

    /// Fill `ebin` from `channels` and `calib_e_fit` (E6).
    pub fn calc_ebins(&mut self) -> Result<(), Error> {
        self.spectrum.ebin = energy_bins(&self.spectrum.channels, &self.calib_e_fit)?;
        Ok(())
    }
}

impl Default for GammaSpectrum {
    fn default() -> Self {
        Self::new()
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn dead_time_subtracts() {
        let mut g = GammaSpectrum::new();
        g.real_time = 209.0;
        g.live_time = 199.0;
        assert_eq!(g.dead_time(), 10.0);
    }

    #[test]
    fn max_channel_empty_is_none() {
        assert_eq!(Spectrum::new("x").max_channel(), None);
    }

    #[test]
    fn calc_ebins_needs_triplet() {
        let mut g = GammaSpectrum::new();
        g.spectrum.channels = vec![0.0, 1.0];
        assert_eq!(g.calc_ebins(), Err(Error::MissingCalibration(0)));
    }

    #[test]
    fn spectrum_len_and_default() {
        let mut s = Spectrum::new("x");
        assert!(s.is_empty());
        assert_eq!(s.len(), 0);
        s.counts = vec![1.0, 2.0];
        s.channels = vec![0.0, 1.0];
        assert!(!s.is_empty());
        assert_eq!(s.len(), 2);
        assert_eq!(s.max_channel(), Some(1.0));
        let g = GammaSpectrum::default();
        assert_eq!(g.spectrum.spec_name, "");
        assert_eq!(g.dead_time(), 0.0);
    }
}
