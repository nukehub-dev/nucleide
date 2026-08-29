import * as React from "react";

export interface LogoProps {
  className?: string;
  size?: number;
  color?: string;
}

export const Logo: React.FC<LogoProps> = ({ className = "", size = 24, color }) => (
  <svg
    xmlns="http://www.w3.org/2000/svg"
    viewBox="0 0 100 100"
    width={size}
    height={size}
    className={className}
    aria-label="Nucleide logo"
    color={color ?? "currentColor"}
  >
    <g fill="currentColor">
      {/* Chart-of-nuclides grid: faint background tiles */}
      <rect x="2.5" y="2.5" width="15" height="15" rx="3" opacity="0.18" />
      <rect x="22.5" y="2.5" width="15" height="15" rx="3" opacity="0.18" />
      <rect x="42.5" y="2.5" width="15" height="15" rx="3" opacity="0.18" />
      <rect x="62.5" y="2.5" width="15" height="15" rx="3" opacity="0.18" />
      <rect x="2.5" y="22.5" width="15" height="15" rx="3" opacity="0.18" />
      <rect x="22.5" y="22.5" width="15" height="15" rx="3" opacity="0.18" />
      <rect x="42.5" y="22.5" width="15" height="15" rx="3" opacity="0.18" />
      <rect x="82.5" y="22.5" width="15" height="15" rx="3" opacity="0.18" />
      <rect x="2.5" y="42.5" width="15" height="15" rx="3" opacity="0.18" />
      <rect x="22.5" y="42.5" width="15" height="15" rx="3" opacity="0.18" />
      <rect x="62.5" y="42.5" width="15" height="15" rx="3" opacity="0.18" />
      <rect x="82.5" y="42.5" width="15" height="15" rx="3" opacity="0.18" />
      <rect x="2.5" y="62.5" width="15" height="15" rx="3" opacity="0.18" />
      <rect x="42.5" y="62.5" width="15" height="15" rx="3" opacity="0.18" />
      <rect x="62.5" y="62.5" width="15" height="15" rx="3" opacity="0.18" />
      <rect x="82.5" y="62.5" width="15" height="15" rx="3" opacity="0.18" />
      <rect x="22.5" y="82.5" width="15" height="15" rx="3" opacity="0.18" />
      <rect x="42.5" y="82.5" width="15" height="15" rx="3" opacity="0.18" />
      <rect x="62.5" y="82.5" width="15" height="15" rx="3" opacity="0.18" />
      <rect x="82.5" y="82.5" width="15" height="15" rx="3" opacity="0.18" />
      {/* Valley-of-stability band (bottom-left to top-right) */}
      <rect x="2.5" y="82.5" width="15" height="15" rx="3" opacity="0.7" />
      <rect x="22.5" y="62.5" width="15" height="15" rx="3" opacity="0.85" />
      <rect x="62.5" y="22.5" width="15" height="15" rx="3" opacity="0.85" />
      <rect x="82.5" y="2.5" width="15" height="15" rx="3" opacity="0.7" />
      {/* Center tile with a punched-out nucleus */}
      <path
        fillRule="evenodd"
        d="M45.5 42.5 h9 a3 3 0 0 1 3 3 v9 a3 3 0 0 1 -3 3 h-9 a3 3 0 0 1 -3 -3 v-9 a3 3 0 0 1 3 -3 z M45.8 50 a4.2 4.2 0 1 1 8.4 0 a4.2 4.2 0 1 1 -8.4 0 z"
      />
    </g>
  </svg>
);
