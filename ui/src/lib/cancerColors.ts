// Consistent color palette for NCI-60 cancer tissue types
// Related biological families share hue families (pink = breast, red = leukemia)
export const CANCER_TYPE_COLORS: Record<string, string> = {
  // Breast family — pink shades
  BREAST:        '#ec4899',
  'MCF7A-REPRO': '#f9a8d4',
  'MCF7D-REPRO': '#be185d',
  // Leukemia family — red shades (NCI-60 repro lines)
  LEUKEMIA:      '#ef4444',
  'K562A-REPRO': '#dc2626',
  'K562B-REPRO': '#f87171',
  // Golub leukemia classes
  AML:           '#ef4444',
  'B-ALL':       '#3b82f6',
  'T-ALL':       '#8b5cf6',
  // Other tissue types — distinct hues
  CNS:           '#8b5cf6',
  COLON:         '#f97316',
  MELANOMA:      '#a3e635',
  NSCLC:         '#3b82f6',
  OVARIAN:       '#f59e0b',
  RENAL:         '#22c55e',
}

export function getCancerColor(cancerType: string, fallbackIndex = 0): string {
  if (CANCER_TYPE_COLORS[cancerType]) return CANCER_TYPE_COLORS[cancerType]
  const FALLBACK = ['#64748b', '#475569', '#334155', '#1e293b']
  return FALLBACK[fallbackIndex % FALLBACK.length]
}
