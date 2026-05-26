import { useEffect, useRef } from "react";
import type { Node } from "reactflow";
import useStore from "../store/useStore";

const PASTE_OFFSET = 50;

type Clipboard = Node[];

export const useCopyPaste = () => {
	const { nodes, setNodes } = useStore();
	const nodesRef = useRef(nodes);
	const clipboard = useRef<Clipboard | null>(null);

	useEffect(() => {
		nodesRef.current = nodes;
	}, [nodes]);

	useEffect(() => {
		const onKeyDown = (e: KeyboardEvent) => {
			if (!(e.ctrlKey || e.metaKey)) return;

			// Do not intercept copy/paste when the user is editing text
			const target = e.target as HTMLElement;
			if (
				target.tagName === "INPUT" ||
				target.tagName === "TEXTAREA" ||
				target.tagName === "SELECT" ||
				target.isContentEditable
			)
				return;

			if (e.key === "c") {
				const selected = nodesRef.current.filter((n) => n.selected);
				if (!selected.length) return;
				clipboard.current = selected;
			}

			if (e.key === "v" && clipboard.current?.length) {
				e.preventDefault();
				const stamp = Date.now().toString(36);
				const newNodes: Node[] = clipboard.current.map((n, i) => {
					// Copy physics/config data; strip identity fields and stale results
					const {
						id: _id,
						label: _label,
						result: _result,
						...physicsData
					} = n.data as Record<string, unknown>;
					const newId = `node_${stamp}_${i}`;
					return {
						...n,
						id: newId,
						position: {
							x: n.position.x + PASTE_OFFSET,
							y: n.position.y + PASTE_OFFSET,
						},
						selected: true,
						data: {
							...physicsData,
							id: newId,
						},
					};
				});

				const current = nodesRef.current;
				setNodes([
					...current.map((n) => ({ ...n, selected: false })),
					...newNodes,
				]);
			}
		};

		document.addEventListener("keydown", onKeyDown);
		return () => document.removeEventListener("keydown", onKeyDown);
	}, [setNodes]);
};
