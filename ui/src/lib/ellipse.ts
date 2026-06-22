/**
 * Compute parametric points along a 2D confidence ellipse from sample coordinates.
 */
export function computeEllipsePoints(
  points: [number, number][],
  nPoints = 60,
  sigmaScale = 1.5,
): [number, number][] {
  if (points.length === 0) return []

  if (points.length === 1) {
    const [x, y] = points[0]
    return Array.from({ length: nPoints + 1 }, (_, i) => {
      const t = (2 * Math.PI * i) / nPoints
      return [x + 0.08 * Math.cos(t), y + 0.08 * Math.sin(t)] as [number, number]
    })
  }

  const n = points.length
  const cx = points.reduce((s, p) => s + p[0], 0) / n
  const cy = points.reduce((s, p) => s + p[1], 0) / n

  let cxx = 0
  let cyy = 0
  let cxy = 0
  for (const [x, y] of points) {
    const dx = x - cx
    const dy = y - cy
    cxx += dx * dx
    cyy += dy * dy
    cxy += dx * dy
  }
  cxx /= n
  cyy /= n
  cxy /= n

  const trace = cxx + cyy
  const det = cxx * cyy - cxy * cxy
  const temp = Math.sqrt(Math.max(0, (trace * trace) / 4 - det))
  const lambda1 = trace / 2 + temp
  const lambda2 = trace / 2 - temp

  const theta = cxy === 0 && cxx >= cyy ? 0 : Math.atan2(lambda1 - cxx, cxy)

  const a = sigmaScale * Math.sqrt(Math.max(lambda1, 0))
  const b = sigmaScale * Math.sqrt(Math.max(lambda2, 0))

  const cosTheta = Math.cos(theta)
  const sinTheta = Math.sin(theta)

  return Array.from({ length: nPoints + 1 }, (_, i) => {
    const t = (2 * Math.PI * i) / nPoints
    const cosT = Math.cos(t)
    const sinT = Math.sin(t)
    const x = cx + a * cosT * cosTheta - b * sinT * sinTheta
    const y = cy + a * cosT * sinTheta + b * sinT * cosTheta
    return [x, y] as [number, number]
  })
}
