/** evaluomeR mark: three found groups huddled on a cobalt tile. Same artwork as public/favicon.svg. */
export function Logo({ size = 32 }: { size?: number }) {
  return (
    <svg className="logo" width={size} height={size} viewBox="0 0 48 48" aria-hidden="true">
      <rect width="48" height="48" rx="13" fill="#1D3FBF" />
      <g stroke="#1D3FBF" strokeWidth="3.2" paintOrder="stroke" strokeLinejoin="round">
        <rect x="8" y="8.5" width="17" height="17" rx="3.2" fill="#FF8FA3" />
        <circle cx="30.5" cy="18.5" r="9.5" fill="#FFD84D" />
        <path d="M23.5 21.5 34.5 40H12.5Z" fill="#7CE0C3" />
      </g>
    </svg>
  )
}
