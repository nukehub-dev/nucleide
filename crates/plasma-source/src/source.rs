//! Source geometry and configuration: point and axisymmetric ring.
//!
//! All lengths are centimetres (the transport-code card convention) and all
//! energies MeV. The v1 geometry set is deliberately small — a point at an
//! arbitrary position and a ring centered on the machine axis (tokamak
//! midplane geometry, uniform azimuth). Eccentric rings and vertical
//! elongation stay out of scope and are rejected loudly
//! ([`Error::NotYetSupported`]) at every entry point that could otherwise
//! guess at them. Partial toroidal sectors live on the parametric plasma
//! config ([`crate::parametric::ToroidalSector`]), not on ring/point
//! sources.

use crate::{Error, FusionReaction, Result, SpectrumSpec};

/// Point source at an arbitrary position.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct PointSource {
    /// x \[cm\].
    pub x_cm: f64,
    /// y \[cm\].
    pub y_cm: f64,
    /// z \[cm\].
    pub z_cm: f64,
}

/// Ring source centered on the machine (z) axis — the tokamak toroidal
/// midplane circle of major radius `radius_cm` at axial height `height_cm`.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct RingSource {
    /// Major radius \[cm\]; finite and strictly positive.
    pub radius_cm: f64,
    /// Axial (vertical) position of the ring plane \[cm\]; finite.
    pub height_cm: f64,
}

/// Source geometry model.
#[derive(Debug, Clone, Copy, PartialEq)]
pub enum SourceModel {
    /// Isotropic point source.
    Point(PointSource),
    /// Axisymmetric ring source (uniform azimuth).
    Ring(RingSource),
}

/// Complete fusion-neutron source configuration.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct PlasmaSourceConfig {
    /// Geometry model.
    pub model: SourceModel,
    /// Fusion reaction (D-D or D-T).
    pub reaction: FusionReaction,
    /// Ion temperature \[keV\] driving the Gaussian broadening; 0 gives the
    /// monoenergetic nominal line.
    pub ion_temperature_kev: f64,
    /// Particle weight carried by the sampler and the emitted cards.
    pub weight: f64,
}

impl PlasmaSourceConfig {
    /// Point source with unit weight.
    pub fn point(
        x_cm: f64,
        y_cm: f64,
        z_cm: f64,
        reaction: FusionReaction,
        ion_temperature_kev: f64,
    ) -> Self {
        Self {
            model: SourceModel::Point(PointSource { x_cm, y_cm, z_cm }),
            reaction,
            ion_temperature_kev,
            weight: 1.0,
        }
    }

    /// Axisymmetric ring source with unit weight.
    pub fn ring(
        radius_cm: f64,
        height_cm: f64,
        reaction: FusionReaction,
        ion_temperature_kev: f64,
    ) -> Self {
        Self {
            model: SourceModel::Ring(RingSource {
                radius_cm,
                height_cm,
            }),
            reaction,
            ion_temperature_kev,
            weight: 1.0,
        }
    }

    /// Override the particle weight (default 1.0).
    pub fn with_weight(mut self, weight: f64) -> Self {
        self.weight = weight;
        self
    }

    /// Validate every field; errors are loud and named.
    pub fn validate(&self) -> Result<()> {
        match self.model {
            SourceModel::Point(p) => {
                for (field, value) in [
                    ("point x", p.x_cm),
                    ("point y", p.y_cm),
                    ("point z", p.z_cm),
                ] {
                    if !value.is_finite() {
                        return Err(Error::NonFinite(field));
                    }
                }
            }
            SourceModel::Ring(r) => {
                if !r.radius_cm.is_finite() {
                    return Err(Error::NonFinite("ring radius"));
                }
                if r.radius_cm <= 0.0 {
                    return Err(Error::NonPositiveRadius(r.radius_cm));
                }
                if !r.height_cm.is_finite() {
                    return Err(Error::NonFinite("ring height"));
                }
            }
        }
        if !self.ion_temperature_kev.is_finite() {
            return Err(Error::NonFinite("ion temperature"));
        }
        if self.ion_temperature_kev < 0.0 {
            return Err(Error::NegativeIonTemperature(self.ion_temperature_kev));
        }
        if !self.weight.is_finite() || self.weight <= 0.0 {
            return Err(Error::NonPositiveWeight(self.weight));
        }
        Ok(())
    }

    /// The temperature-broadened spectrum of this configuration.
    pub fn spectrum(&self) -> Result<SpectrumSpec> {
        self.validate()?;
        SpectrumSpec::from_reaction(self.reaction, self.ion_temperature_kev)
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn validation_accepts_conventional_sources() {
        PlasmaSourceConfig::point(1.0, 2.0, 3.0, FusionReaction::Dt, 20.0)
            .validate()
            .unwrap();
        PlasmaSourceConfig::ring(300.0, 0.0, FusionReaction::Dd, 0.0)
            .validate()
            .unwrap();
        PlasmaSourceConfig::ring(300.0, -40.0, FusionReaction::Dt, 10.0)
            .with_weight(0.5)
            .validate()
            .unwrap();
    }

    #[test]
    fn invalid_fields_are_loud() {
        assert_eq!(
            PlasmaSourceConfig::point(f64::NAN, 0.0, 0.0, FusionReaction::Dt, 0.0).validate(),
            Err(Error::NonFinite("point x"))
        );
        assert_eq!(
            PlasmaSourceConfig::ring(0.0, 0.0, FusionReaction::Dt, 0.0).validate(),
            Err(Error::NonPositiveRadius(0.0))
        );
        assert_eq!(
            PlasmaSourceConfig::ring(-5.0, 0.0, FusionReaction::Dt, 0.0).validate(),
            Err(Error::NonPositiveRadius(-5.0))
        );
        assert_eq!(
            PlasmaSourceConfig::ring(f64::INFINITY, 0.0, FusionReaction::Dt, 0.0).validate(),
            Err(Error::NonFinite("ring radius"))
        );
        assert_eq!(
            PlasmaSourceConfig::ring(300.0, f64::NAN, FusionReaction::Dt, 0.0).validate(),
            Err(Error::NonFinite("ring height"))
        );
        assert_eq!(
            PlasmaSourceConfig::ring(300.0, 0.0, FusionReaction::Dt, -1.0).validate(),
            Err(Error::NegativeIonTemperature(-1.0))
        );
        assert_eq!(
            PlasmaSourceConfig::ring(300.0, 0.0, FusionReaction::Dt, 0.0)
                .with_weight(0.0)
                .validate(),
            Err(Error::NonPositiveWeight(0.0))
        );
    }

    #[test]
    fn spectrum_builds_through_validation() {
        let config = PlasmaSourceConfig::ring(300.0, 0.0, FusionReaction::Dt, 20.0);
        let spec = config.spectrum().unwrap();
        assert!(!spec.is_mono());
        let bad = PlasmaSourceConfig::ring(300.0, 0.0, FusionReaction::Dt, -2.0);
        assert_eq!(bad.spectrum(), Err(Error::NegativeIonTemperature(-2.0)));
    }
}
