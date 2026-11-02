import type { DemoData } from '@/types/demo'

// Cluster colours, readable on cobalt. Each group also gets a shape so colour is never the only cue.
export const COLORS = ['#FFD84D', '#FF8FA3', '#7CE0C3', '#C9B6FF', '#FFB067', '#8FD3FF', '#F2F2F2', '#B5E36B', '#FF6B6B', '#9AA6FF', '#E0A3FF', '#66D9E8']
export const HAZE = '#C3CEF8'

const SHAPES = [
  (r: number) => `M${-r},0a${r},${r} 0 1,0 ${2 * r},0a${r},${r} 0 1,0 ${-2 * r},0`,
  (r: number) => `M${-r * 0.9},${-r * 0.9}h${r * 1.8}v${r * 1.8}h${-r * 1.8}z`,
  (r: number) => `M0,${-r * 1.1}L${r * 1.05},${r * 0.8}H${-r * 1.05}z`,
  (r: number) => `M0,${-r * 1.2}L${r * 1.1},0L0,${r * 1.2}L${-r * 1.1},0z`,
  (r: number) => `M0,${r * 1.1}L${r * 1.05},${-r * 0.8}H${-r * 1.05}z`,
  (r: number) => `M${-r * 0.4},${-r}h${r * 0.8}v${r * 0.6}h${r * 0.6}v${r * 0.8}h${-r * 0.6}v${r * 0.6}h${-r * 0.8}v${-r * 0.6}h${-r * 0.6}v${-r * 0.8}h${r * 0.6}z`,
]
export const shape = (group: number, r: number) => SHAPES[group % SHAPES.length](r)
export const color = (group: number) => COLORS[group % COLORS.length]

export interface Dataset {
  key: string
  label: string
  short: string
  noun: string
  data: DemoData
}

export const kValues = (d: DemoData) => Object.keys(d.kSummary).map(Number).sort((a, b) => a - b)

/** Maps PC1/PC2 into a [x0,x1]×[y0,y1] box (y up). */
export function project(d: DemoData, box: [number, number, number, number]) {
  const xs = d.pca.coords.map((c) => c.pc1)
  const ys = d.pca.coords.map((c) => c.pc2)
  const [mx, Mx, my, My] = [Math.min(...xs), Math.max(...xs), Math.min(...ys), Math.max(...ys)]
  const [x0, x1, y0, y1] = box
  return d.pca.coords.map((c) => ({
    x: x0 + ((c.pc1 - mx) / (Mx - mx)) * (x1 - x0),
    y: y1 - ((c.pc2 - my) / (My - my)) * (y1 - y0),
  }))
}
