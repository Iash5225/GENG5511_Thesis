def extract_datasets_with_meta(data):
    """
    Reuses extract_datasets() for numeric arrays, then builds meta dict
    (Author, Year) per property key. No duplication of numeric work.
    """
    datasets = extract_datasets(data)

    def meta_block(df, n):
        if df is None:
            return None
        authors = df.get("Author", pd.Series(["Unknown"] * n)).to_numpy()
        years = df.get("Year", pd.Series([0] * n))
        # Ensure year is integer (handle non-numeric gracefully)
        years = pd.to_numeric(years, errors="coerce").fillna(
            0).astype(int).to_numpy()
        return {"Author": authors, "Year": years}

    meta = {
        "Vm_sub":       meta_block(data.get("cell_volume_sub"),   len(datasets[0])),
        "Vm_melt":      meta_block(data.get("cell_volume_melt"),  len(datasets[3])),
        "Vm_highp":     meta_block(data.get("cell_volume_highp"), len(datasets[6])) if len(datasets[6]) else None,
        "cp_sub":       meta_block(data.get("heat_capacity"),     len(datasets[9])),
        "alpha_sub":    meta_block(data.get("thermal_coeff"),     len(datasets[12])),
        "BetaT_sub":    meta_block(data.get("bulk_t"),            len(datasets[15])),
        "BetaS_sub":    meta_block(data.get("bulk_s"),            len(datasets[18])),
        "sublimation":  meta_block(data.get("sublimation"),       len(datasets[21])),
        "melting":      meta_block(data.get("melting"),           len(datasets[26])),
        "H_solid_sub":  meta_block(data.get("heatsub"),           len(datasets[31])),
        "H_solid_melt": meta_block(data.get("fusion"),            len(datasets[34])),
    }
    return datasets, meta
