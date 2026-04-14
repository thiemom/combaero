feat(gui): integrate surface enhancements and verify analytical Jacobians

- Implement fixed/verified C++ analytical Jacobians for ribbed, dimpled, pin-fin, and impingement models.
- Extend gui backend schemas and graph builder to support modular surface enhancements.
- Create SurfaceEnhancementInspector component for dynamic parameter configuration in React frontend.
- Enrich node and element diagnostics to surface Nu, htc, and dP metrics for physical verification.
- Pass all 1050+ backend tests including newly added GUI serialization logic.
