use semver::Version;

pub const VALID_VERSIONS_PACKAGE: &[Version] = &[Version::new(0, 1, 0), Version::new(0, 2, 0)];
pub const VALID_VERSIONS_PITCH: &[Version] = &[Version::new(0, 1, 0)];
pub const VALID_VERSIONS_PROPAGATION: &[Version] = &[Version::new(0, 1, 0), Version::new(0, 2, 0)];
pub const VALID_VERSIONS_QA: &[Version] = &[Version::new(0, 1, 0)];
pub const VALID_VERSIONS_REVIVE: &[Version] = &[Version::new(0, 1, 0)];
pub const VALID_VERSIONS_UNIFORMITY: &[Version] = &[
    Version::new(0, 1, 0),
    Version::new(0, 2, 0),
    Version::new(0, 2, 1),
];
pub const VALID_VERSIONS_UNIFORMITY_MISSING_TEMPERATURE: &[Version] = &[Version::new(0, 1, 0)];

pub const VALID_VERSIONS_PLATING: &[Version] = &[Version::new(0, 1, 0)];
pub const VALID_VERSIONS_THERMAL: &[Version] = &[Version::new(0, 1, 0)];
