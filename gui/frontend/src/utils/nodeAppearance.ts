/**
 * One icon and one colour per node type (#321 item 5).
 *
 * The palette (Sidebar) and the canvas cards (components/nodes/*) used to pick
 * their own, independently: 7 of 17 types drifted apart and three pairs shared
 * a glyph. Both now read this table, so they cannot diverge again. A canvas
 * card may still choose a variant from data (an area change draws Minimize2
 * when it contracts) but takes its colour from here.
 */

import {
	Activity,
	ArrowRight,
	BetweenHorizontalStart,
	CircleDot,
	Columns4,
	Crosshair,
	Database,
	Flame,
	Gauge,
	Grip,
	Heater,
	Link,
	LogIn,
	type LucideIcon,
	Maximize2,
	ShowerHead,
	Split,
	Square,
	Tally2,
	Tornado,
	Waves,
	Wind,
} from "lucide-react";

export interface NodeAppearance {
	icon: LucideIcon;
	/** Tailwind text colour of the glyph. */
	colour: string;
	/** Extra classes some glyphs need (a filled square). */
	extra?: string;
}

export const NODE_APPEARANCE = {
	mass_boundary: { icon: LogIn, colour: "text-orange-600" },
	pressure_boundary: { icon: Gauge, colour: "text-blue-600" },
	plenum: { icon: Database, colour: "text-amber-600" },
	wall: { icon: Square, colour: "text-stone-600", extra: "fill-stone-600" },
	combustor: { icon: Flame, colour: "text-red-600" },
	momentum_chamber: { icon: CircleDot, colour: "text-indigo-600" },
	channel: { icon: ArrowRight, colour: "text-stone-600" },
	orifice: { icon: Tally2, colour: "text-orange-500" },
	impingement_array: { icon: ShowerHead, colour: "text-orange-500" },
	effusion_plate: { icon: Grip, colour: "text-orange-500" },
	film_cooling: { icon: Waves, colour: "text-orange-500" },
	slot_film: { icon: BetweenHorizontalStart, colour: "text-orange-500" },
	effusion_liner: { icon: Columns4, colour: "text-orange-500" },
	area_change: { icon: Maximize2, colour: "text-blue-500" },
	discrete_loss: { icon: Activity, colour: "text-purple-500" },
	lossless_connection: { icon: Link, colour: "text-slate-500" },
	mpce_tee: { icon: Split, colour: "text-pink-600" },
	ejector: { icon: Wind, colour: "text-sky-500" },
	vortex: { icon: Tornado, colour: "text-violet-500" },
	probe: { icon: Crosshair, colour: "text-blue-500" },
	heat_source: { icon: Heater, colour: "text-rose-600" },
} satisfies Record<string, NodeAppearance>;

export type AppearanceType = keyof typeof NODE_APPEARANCE;

/** The glyph's className: its colour plus any extra classes. */
export function iconClass(type: AppearanceType, more = ""): string {
	const a: NodeAppearance = NODE_APPEARANCE[type];
	return [a.colour, a.extra, more].filter(Boolean).join(" ");
}
