interface GridMargins {
  left: number
  right: number
  top: number
}

function animatedLineSegment(kValue: number | string, yMin: number, yMax: number, lineColor: string) {
  return [
    {
      coord: [String(kValue), yMin],
      symbol: 'none' as const,
      lineStyle: { color: lineColor, type: 'dashed' as const, width: 1.5 },
      label: { show: false },
    },
    {
      coord: [String(kValue), yMax],
      symbol: 'none' as const,
      label: { show: false },
    },
  ]
}

function columnXPercent(index: number, count: number, grid: GridMargins) {
  const plotWidth = 100 - grid.left - grid.right
  return grid.left + ((index + 0.5) / count) * plotWidth
}

function graphicTextElement(
  categoryKey: string,
  categories: string[],
  grid: GridMargins,
  text: string,
  color: string,
  bold = false,
) {
  const index = categories.indexOf(categoryKey)
  if (index < 0) return null

  const top = `${Math.max(0, grid.top - 2)}%`
  const left = `${columnXPercent(index, categories.length, grid)}%`

  return {
    type: 'text' as const,
    left,
    top,
    z: 100,
    style: {
      text,
      fill: color,
      fontSize: 10,
      fontWeight: bold ? ('bold' as const) : ('normal' as const),
      textAlign: 'center' as const,
      backgroundColor: 'rgba(15,23,42,0.7)',
      padding: [2, 6] as [number, number],
      borderRadius: 4,
    },
  }
}

export function buildAnimatedLines(
  selectedK: number,
  optimalK: number,
  yMin: number,
  yMax: number,
) {
  return {
    markLine: {
      silent: true,
      symbol: 'none',
      animation: true,
      animationDuration: 500,
      animationEasing: 'cubicOut' as const,
      data: [
        animatedLineSegment(optimalK, yMin, yMax, '#22c55e'),
        ...(selectedK !== optimalK
          ? [animatedLineSegment(selectedK, yMin, yMax, '#cbd5e1')]
          : []),
      ],
    },
  }
}

export function buildGraphicKLabels(
  selectedK: number,
  optimalK: number,
  categories: string[],
  grid: GridMargins,
) {
  const elements = [
    graphicTextElement(
      String(optimalK),
      categories,
      grid,
      `✦ Optimal k = ${optimalK}`,
      '#22c55e',
      true,
    ),
    ...(selectedK !== optimalK
      ? [
          graphicTextElement(
            String(selectedK),
            categories,
            grid,
            `k = ${selectedK}`,
            '#cbd5e1',
          ),
        ]
      : []),
  ].filter((el): el is NonNullable<typeof el> => el !== null)

  return { graphic: elements }
}

export function buildOptimalAnimatedLine(kValue: number | string, yMin: number, yMax: number) {
  return {
    markLine: {
      silent: true,
      symbol: 'none',
      animation: true,
      animationDuration: 500,
      animationEasing: 'cubicOut' as const,
      data: [animatedLineSegment(kValue, yMin, yMax, '#22c55e')],
    },
  }
}

export function buildGraphicOptimalLabel(
  kValue: number | string,
  categories: string[],
  grid: GridMargins,
  labelFormatter: string,
) {
  const element = graphicTextElement(String(kValue), categories, grid, labelFormatter, '#22c55e', true)
  return { graphic: element ? [element] : [] }
}
