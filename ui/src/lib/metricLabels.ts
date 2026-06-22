/** Human-readable label for evaluomeR metric / index names */
export function formatMetricLabel(name: string): string {
  if (name === 'all_metrics') return 'All metrics'
  return name
}
