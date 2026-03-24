use anyhow::{bail, Result};
use std::collections::HashMap;

/// A clone in the tumour's clonal architecture.
#[derive(Debug, Clone)]
pub struct Clone {
    pub id: String,
    pub ccf: f64,
    pub parent: Option<String>,
}

/// Represents the clonal tree of a tumour.
///
/// The root clone has CCF = 1.0 (all tumour cells carry its mutations).
/// Subclones have CCF <= parent CCF.
#[derive(Debug, Clone)]
pub struct ClonalTree {
    clones: Vec<Clone>,
    children: HashMap<String, Vec<String>>,
}

impl ClonalTree {
    /// Build and validate a clonal tree from a list of clones.
    pub fn new(clones: Vec<Clone>) -> Result<Self> {
        if clones.is_empty() {
            bail!("clonal tree must have at least one clone");
        }

        let ccf_map: HashMap<&str, f64> = clones.iter().map(|c| (c.id.as_str(), c.ccf)).collect();

        // Validate CCFs are in [0, 1]
        for clone in &clones {
            if !(0.0..=1.0).contains(&clone.ccf) {
                bail!(
                    "clone {} CCF {} must be between 0.0 and 1.0",
                    clone.id,
                    clone.ccf
                );
            }
        }

        // Validate parent references exist and child CCF <= parent CCF
        let mut children: HashMap<String, Vec<String>> = HashMap::new();
        for clone in &clones {
            if let Some(parent_id) = &clone.parent {
                let parent_ccf = ccf_map.get(parent_id.as_str()).ok_or_else(|| {
                    anyhow::anyhow!("clone {} references unknown parent {}", clone.id, parent_id)
                })?;
                if clone.ccf > *parent_ccf {
                    bail!(
                        "clone {} CCF ({}) exceeds parent {} CCF ({})",
                        clone.id,
                        clone.ccf,
                        parent_id,
                        parent_ccf
                    );
                }
                children
                    .entry(parent_id.clone())
                    .or_default()
                    .push(clone.id.clone());
            }
        }

        // Validate exactly one root (no parent)
        let roots: Vec<_> = clones.iter().filter(|c| c.parent.is_none()).collect();
        if roots.len() != 1 {
            bail!(
                "clonal tree must have exactly one root clone, found {}",
                roots.len()
            );
        }

        Ok(Self { clones, children })
    }

    /// Get all clones in the tree.
    pub fn clones(&self) -> &[Clone] {
        &self.clones
    }

    /// Get the root clone.
    #[allow(dead_code)]
    pub fn root(&self) -> &Clone {
        self.clones.iter().find(|c| c.parent.is_none()).unwrap()
    }

    /// Get children of a clone.
    pub fn children_of(&self, clone_id: &str) -> &[String] {
        self.children
            .get(clone_id)
            .map(|v| v.as_slice())
            .unwrap_or(&[])
    }

    /// Get a clone by ID.
    pub fn get(&self, clone_id: &str) -> Option<&Clone> {
        self.clones.iter().find(|c| c.id == clone_id)
    }

    /// Get all ancestor clone IDs for a given clone (including itself).
    /// Mutations are inherited from all ancestors.
    #[allow(dead_code)]
    pub fn ancestors(&self, clone_id: &str) -> Vec<String> {
        let mut result = Vec::new();
        let mut current = Some(clone_id.to_string());
        while let Some(id) = current {
            result.push(id.clone());
            current = self.get(&id).and_then(|c| c.parent.clone());
        }
        result.reverse();
        result
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn simple_tree() -> ClonalTree {
        ClonalTree::new(vec![
            Clone {
                id: "root".into(),
                ccf: 1.0,
                parent: None,
            },
            Clone {
                id: "sub_a".into(),
                ccf: 0.4,
                parent: Some("root".into()),
            },
            Clone {
                id: "sub_b".into(),
                ccf: 0.2,
                parent: Some("root".into()),
            },
            Clone {
                id: "sub_a1".into(),
                ccf: 0.1,
                parent: Some("sub_a".into()),
            },
        ])
        .unwrap()
    }

    #[test]
    fn test_valid_tree() {
        let tree = simple_tree();
        assert_eq!(tree.clones().len(), 4);
        assert_eq!(tree.root().id, "root");
    }

    #[test]
    fn test_children() {
        let tree = simple_tree();
        let children = tree.children_of("root");
        assert_eq!(children.len(), 2);
        assert!(children.contains(&"sub_a".to_string()));
        assert!(children.contains(&"sub_b".to_string()));
    }

    #[test]
    fn test_ancestors() {
        let tree = simple_tree();
        let ancestors = tree.ancestors("sub_a1");
        assert_eq!(ancestors, vec!["root", "sub_a", "sub_a1"]);
    }

    #[test]
    fn test_root_ancestors() {
        let tree = simple_tree();
        let ancestors = tree.ancestors("root");
        assert_eq!(ancestors, vec!["root"]);
    }

    #[test]
    fn test_reject_child_ccf_exceeds_parent() {
        let result = ClonalTree::new(vec![
            Clone {
                id: "root".into(),
                ccf: 0.5,
                parent: None,
            },
            Clone {
                id: "sub".into(),
                ccf: 0.8,
                parent: Some("root".into()),
            },
        ]);
        assert!(result.is_err());
    }

    #[test]
    fn test_reject_unknown_parent() {
        let result = ClonalTree::new(vec![
            Clone {
                id: "root".into(),
                ccf: 1.0,
                parent: None,
            },
            Clone {
                id: "sub".into(),
                ccf: 0.5,
                parent: Some("nonexistent".into()),
            },
        ]);
        assert!(result.is_err());
    }

    #[test]
    fn test_reject_no_root() {
        let result = ClonalTree::new(vec![
            Clone {
                id: "a".into(),
                ccf: 0.5,
                parent: Some("b".into()),
            },
            Clone {
                id: "b".into(),
                ccf: 0.8,
                parent: Some("a".into()),
            },
        ]);
        assert!(result.is_err());
    }

    #[test]
    fn test_reject_multiple_roots() {
        let result = ClonalTree::new(vec![
            Clone {
                id: "a".into(),
                ccf: 1.0,
                parent: None,
            },
            Clone {
                id: "b".into(),
                ccf: 1.0,
                parent: None,
            },
        ]);
        assert!(result.is_err());
    }

    #[test]
    fn test_reject_empty() {
        let result = ClonalTree::new(vec![]);
        assert!(result.is_err());
    }

    #[test]
    fn test_ccf_out_of_range() {
        let result = ClonalTree::new(vec![Clone {
            id: "root".into(),
            ccf: 1.5,
            parent: None,
        }]);
        assert!(result.is_err());
    }

    #[test]
    fn test_single_clone() {
        let tree = ClonalTree::new(vec![Clone {
            id: "root".into(),
            ccf: 1.0,
            parent: None,
        }])
        .unwrap();
        assert_eq!(tree.clones().len(), 1);
        assert!(tree.children_of("root").is_empty());
    }
}
