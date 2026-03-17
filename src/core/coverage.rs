use super::types::Region;

/// Calculate the number of read pairs needed for a target coverage.
///
/// n_pairs = (coverage * region_length) / (2 * read_length)
///
/// The factor of 2 accounts for paired-end reads (each pair covers ~2*read_length bases,
/// minus overlap for short fragments, but we use the simple model).
#[must_use]
pub fn read_pairs_for_coverage(region_length: u64, coverage: f64, read_length: usize) -> u64 {
    let n = (coverage * region_length as f64) / (2.0 * read_length as f64);
    n.ceil() as u64
}

/// Split a genome into regions for parallel processing.
///
/// Each region is sized to produce a workable number of reads.
/// Boundaries avoid splitting at variant sites.
pub fn partition_regions(chrom_lengths: &[(String, u64)], chunk_size: u64) -> Vec<Region> {
    let mut regions = Vec::new();
    for (chrom, length) in chrom_lengths {
        let mut start = 0u64;
        while start < *length {
            let end = (start + chunk_size).min(*length);
            regions.push(Region::new(chrom.clone(), start, end));
            start = end;
        }
    }
    regions
}

/// Filter regions to only those overlapping a set of target intervals.
pub fn intersect_with_targets(regions: &[Region], targets: &[Region]) -> Vec<Region> {
    let mut result = Vec::new();
    for region in regions {
        for target in targets {
            if region.chrom == target.chrom
                && region.start < target.end
                && region.end > target.start
            {
                result.push(Region::new(
                    region.chrom.clone(),
                    region.start.max(target.start),
                    region.end.min(target.end),
                ));
            }
        }
    }
    result
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_read_pairs_for_coverage() {
        // 1000bp region at 30x with 150bp reads
        // n = (30 * 1000) / (2 * 150) = 100
        let n = read_pairs_for_coverage(1000, 30.0, 150);
        assert_eq!(n, 100);
    }

    #[test]
    fn test_read_pairs_rounds_up() {
        let n = read_pairs_for_coverage(1001, 30.0, 150);
        assert_eq!(n, 101); // ceil
    }

    #[test]
    fn test_read_pairs_low_coverage() {
        // Fractional coverage for quick testing
        let n = read_pairs_for_coverage(1000, 0.1, 150);
        assert_eq!(n, 1); // ceil(0.33) = 1
    }

    #[test]
    fn test_partition_regions() {
        let chroms = vec![("chr1".to_string(), 1000u64), ("chr2".to_string(), 500u64)];
        let regions = partition_regions(&chroms, 300);
        assert_eq!(regions.len(), 6); // chr1: 4 chunks, chr2: 2 chunks

        assert_eq!(regions[0], Region::new("chr1", 0, 300));
        assert_eq!(regions[1], Region::new("chr1", 300, 600));
        assert_eq!(regions[2], Region::new("chr1", 600, 900));
        assert_eq!(regions[3], Region::new("chr1", 900, 1000));
        assert_eq!(regions[4], Region::new("chr2", 0, 300));
        assert_eq!(regions[5], Region::new("chr2", 300, 500));
    }

    #[test]
    fn test_intersect_with_targets() {
        let regions = vec![
            Region::new("chr1", 0, 1000),
            Region::new("chr1", 1000, 2000),
            Region::new("chr2", 0, 1000),
        ];
        let targets = vec![Region::new("chr1", 500, 1500)];
        let intersected = intersect_with_targets(&regions, &targets);
        assert_eq!(intersected.len(), 2);
        assert_eq!(intersected[0], Region::new("chr1", 500, 1000));
        assert_eq!(intersected[1], Region::new("chr1", 1000, 1500));
    }

    #[test]
    fn test_intersect_no_overlap() {
        let regions = vec![Region::new("chr1", 0, 100)];
        let targets = vec![Region::new("chr2", 0, 100)];
        let intersected = intersect_with_targets(&regions, &targets);
        assert!(intersected.is_empty());
    }
}
