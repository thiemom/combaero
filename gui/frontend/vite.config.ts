import react from "@vitejs/plugin-react";
import { defineConfig } from "vite";

// https://vite.dev/config/
export default defineConfig({
	plugins: [react()],
	server: {
		proxy: {
			"/solve": "http://localhost:8000",
			"/export": "http://localhost:8000",
			"/metadata": "http://localhost:8000",
			"/solver": "http://localhost:8000",
			"/health": "http://localhost:8000",
		},
	},
});
